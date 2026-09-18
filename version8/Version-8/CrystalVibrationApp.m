classdef CrystalVibrationApp < handle
    %CRYSTALVIBRATIONAPP GUI: load a CIF/txt crystal, build the spring network, solve the modes and animate them.
    %   Build & Solve: autoRigidShells/generateLatticeGeneral -> buildDynamicalMatrix -> jacobiEigenSolver.
    %   Views window: 3D lattice (hover an atom for its element, via atomHover + getElementDetails) and frequency spectrum.
    %   Shear Analysis window: slider (applyShearForce), Relax button (relaxShearLattice),
    %   and the V7 buttons that run phononShearSuite and open plotShearSoftening / plotShearBonds figures.

    properties (Access = private)
        UIFigure
        LoadFileButton
        StatusLabel
        NEditField
        kEditField
        cutoffDropDown
        BuildButton
        ModeDropDown
        AmplitudeSlider
        SpeedSlider
        PlayButton
        PauseButton
        StopButton
        LaunchShearButton

        PlotFig
        Panel3D
        PanelSpec
        PlotAx
        SpecAx
        View3DButton
        ViewSpecButton
        SpecHighlightHandle

        ShearFig
        ShearGrid
        ShearControlPanel
        ShearLabel
        ShearSlider
        RelaxButton
        RelaxStatusLabel
        UIAxes3D_NoForce
        UIAxes3D_WithForce
        UIAxesMode_NoForce
        UIAxesMode_WithForce

        PhononPlotsButton
        BondPlotsButton
        DispersionButton
        GammaMaxEdit
        AlphaEdit
        AlphaLabel
        AnalysisStatusLabel
        AnalysisLV
        AnalysisBasis
        AnalysisCache

        latticeVectors
        basisFrac
        basisMasses
        useCustomFile = false

        pos
        bonds
        masses
        V
        omega
        isBuilt = false
        BaseStiffnessMatrix
        BaseSpringConstant

        modeShape
        currentOmega = 0
        animTime = 0
        animTimer
        atomHandle
        bondHandles
    end

    properties (Constant, Access = private)
        AnimFrameRateHz = 30
    end

    methods (Access = public)
        function app = CrystalVibrationApp()
            app.createComponents();
        end

        function delete(app)
            app.stopTimerIfRunning();
            if ~isempty(app.PlotFig) && isvalid(app.PlotFig)
                delete(app.PlotFig);
            end
            if ~isempty(app.ShearFig) && isvalid(app.ShearFig)
                delete(app.ShearFig);
            end
            if isvalid(app.UIFigure)
                delete(app.UIFigure);
            end
        end
    end

    methods (Access = private)

        function createComponents(app)
            app.UIFigure = uifigure('Name', 'Crystal Vibration Explorer - Controls', ...
                'Position', [100 100 320 650], 'Color', 'w');
            app.UIFigure.CloseRequestFcn = @(~, ~) app.onClose();

            uilabel(app.UIFigure, 'Text', 'Crystal Vibration Explorer', ...
                'FontSize', 16, 'FontWeight', 'bold', 'Position', [20 610 280 25]);

            app.LoadFileButton = uibutton(app.UIFigure, 'push', ...
                'Text', 'Load Crystal File...', 'Position', [20 575 260 30], ...
                'ButtonPushedFcn', @(~, ~) app.LoadFileButtonPushed());

            app.StatusLabel = uilabel(app.UIFigure, 'Text', 'No file loaded (using default bcc structure)', ...
                'Position', [20 550 280 20], 'FontColor', [0.4 0.4 0.4]);

            uilabel(app.UIFigure, 'Text', 'Supercell size (fixed): 1 x 1 x 1 cell', 'Position', [20 515 260 20]);
            app.NEditField = uieditfield(app.UIFigure, 'numeric', ...
                'Value', 1, 'Limits', [1 1], 'Editable', 'off', 'Position', [20 490 260 25]);

            uilabel(app.UIFigure, 'Text', 'Spring constant k:', 'Position', [20 460 260 20]);
            % Spring constant k is a reduced-unit default (1.0): use a force constant fitted to the real material.
            app.kEditField = uieditfield(app.UIFigure, 'numeric', ...
                'Value', 1.0, 'Limits', [0.001 Inf], 'Position', [20 435 260 25]);

            uilabel(app.UIFigure, 'Text', 'Bond network:', 'Position', [20 400 260 20]);
            app.cutoffDropDown = uidropdown(app.UIFigure, ...
                'Items', {'Nearest neighbors only', '1st + 2nd neighbors (fixed)', 'Auto (grows shells until rigid)'}, ...
                'Value', 'Auto (grows shells until rigid)', ...
                'Position', [20 375 260 25]);

            app.BuildButton = uibutton(app.UIFigure, 'push', ...
                'Text', 'Build && Solve', 'Position', [20 330 260 35], ...
                'FontWeight', 'bold', 'BackgroundColor', [0.20 0.45 0.80], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.BuildButtonPushed());

            uilabel(app.UIFigure, 'Text', 'Vibrational mode:', 'Position', [20 295 260 20]);
            app.ModeDropDown = uidropdown(app.UIFigure, ...
                'Items', {'(build the model first)'}, 'Enable', 'off', ...
                'Position', [20 270 260 25], ...
                'ValueChangedFcn', @(~, ~) app.ModeDropDownValueChanged());

            uilabel(app.UIFigure, 'Text', 'Amplitude:', 'Position', [20 235 260 20]);
            app.AmplitudeSlider = uislider(app.UIFigure, ...
                'Limits', [0 1], 'Value', 0.25, 'Position', [25 220 250 3]);

            uilabel(app.UIFigure, 'Text', 'Playback speed (cycles/sec):', 'Position', [20 175 260 20]);
            app.SpeedSlider = uislider(app.UIFigure, ...
                'Limits', [0.1 10], 'Value', 2, 'Position', [25 160 250 3], ...
                'ValueChangedFcn', @(~, ~) app.SpeedSliderValueChanged());

            app.PlayButton = uibutton(app.UIFigure, 'push', 'Text', 'Play', ...
                'Position', [20 105 80 32], 'Enable', 'off', ...
                'BackgroundColor', [0.30 0.70 0.35], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.PlayButtonPushed());
            app.PauseButton = uibutton(app.UIFigure, 'push', 'Text', 'Pause', ...
                'Position', [110 105 80 32], 'Enable', 'off', ...
                'ButtonPushedFcn', @(~, ~) app.PauseButtonPushed());
            app.StopButton = uibutton(app.UIFigure, 'push', 'Text', 'Stop', ...
                'Position', [200 105 80 32], 'Enable', 'off', ...
                'ButtonPushedFcn', @(~, ~) app.StopButtonPushed());

            app.LaunchShearButton = uibutton(app.UIFigure, 'push', ...
                'Text', 'Open Shear Analysis...', 'Position', [20 55 260 35], ...
                'Enable', 'off', 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.4 0.2], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.ensureShearFigure());

            uilabel(app.UIFigure, 'Text', 'The 3D view and spectrum open in a separate window.', ...
                'Position', [20 15 280 30], 'FontColor', [0.4 0.4 0.4]);

            app.ensurePlotFigure();
        end

        function ensurePlotFigure(app)
            if ~isempty(app.PlotFig) && isvalid(app.PlotFig)
                return;
            end

            app.PlotFig = figure('Name', 'Crystal Vibration Explorer - Views', ...
                'NumberTitle', 'off', 'Position', [440 80 680 700], 'Color', 'w');

            app.View3DButton = uicontrol(app.PlotFig, 'Style', 'pushbutton', ...
                'String', '3D View', 'Position', [10 654 200 34], ...
                'FontWeight', 'bold', 'Callback', @(~, ~) app.switchToView('3D'));
            app.ViewSpecButton = uicontrol(app.PlotFig, 'Style', 'pushbutton', ...
                'String', 'Frequency Spectrum', 'Position', [220 654 200 34], ...
                'FontWeight', 'bold', 'Callback', @(~, ~) app.switchToView('Spec'));

            app.Panel3D = uipanel(app.PlotFig, 'Units', 'pixels', ...
                'Position', [10 10 660 630], 'BorderType', 'none', 'BackgroundColor', 'w');
            app.PlotAx = axes('Parent', app.Panel3D);
            title(app.PlotAx, 'Load a crystal or click "Build & Solve" to begin');
            axis(app.PlotAx, 'equal');
            grid(app.PlotAx, 'on');
            view(app.PlotAx, 3);

            app.PanelSpec = uipanel(app.PlotFig, 'Units', 'pixels', ...
                'Position', [10 10 660 630], 'BorderType', 'none', 'BackgroundColor', 'w');
            app.SpecAx = axes('Parent', app.PanelSpec);
            title(app.SpecAx, 'Click "Build & Solve" to compute the spectrum');
            xlabel(app.SpecAx, 'Mode index');
            ylabel(app.SpecAx, '\omega (angular frequency)');
            grid(app.SpecAx, 'on');

            app.switchToView('3D');
        end

        function ensureShearFigure(app)
            if ~app.isBuilt
                app.showAlert('Please Build & Solve the lattice before running shear analysis.', 'Not Built');
                return;
            end

            if ~isempty(app.ShearFig) && isvalid(app.ShearFig)
                figure(app.ShearFig);
                return;
            end

            app.ShearFig = uifigure('Name', 'Crystal Vibration Explorer - Shear Analysis');
            app.ShearFig.Position = [100 100 1400 850];

            app.ShearGrid = uigridlayout(app.ShearFig);
            app.ShearGrid.ColumnWidth = {350, '1x', '1x'};
            app.ShearGrid.RowHeight = {'1x', '1x'};

            app.ShearControlPanel = uipanel(app.ShearGrid);
            app.ShearControlPanel.Title = 'Shear Control Panel';
            app.ShearControlPanel.Layout.Row = [1 2];
            app.ShearControlPanel.Layout.Column = 1;

            app.ShearLabel = uilabel(app.ShearControlPanel);
            app.ShearLabel.Position = [20 750 220 22];
            app.ShearLabel.Text = 'Engineering Shear Strain (gamma):';

            app.ShearSlider = uislider(app.ShearControlPanel);
            app.ShearSlider.Limits = [0 0.5];
            app.ShearSlider.Position = [20 720 300 3];
            app.ShearSlider.ValueChangedFcn = @(src, event) app.updateShearAnalysis(event);

            app.RelaxButton = uibutton(app.ShearControlPanel, 'push', ...
                'Text', 'Relax (Gauss-Jordan, Exp.5)', 'Position', [20 670 300 32], ...
                'FontWeight', 'bold', 'BackgroundColor', [0.20 0.45 0.80], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.RelaxButtonPushed());

            app.RelaxStatusLabel = uilabel(app.ShearControlPanel);
            app.RelaxStatusLabel.Position = [20 630 310 34];
            app.RelaxStatusLabel.Text = 'Relax to solve for true equilibrium of interior atoms.';
            app.RelaxStatusLabel.FontColor = [0.4 0.4 0.4];
            app.RelaxStatusLabel.WordWrap = 'on';

            hdr = uilabel(app.ShearControlPanel, 'Text', 'Phonon softening & bond breaking (V7)', ...
                'Position', [20 590 310 22], 'FontWeight', 'bold');
            uilabel(app.ShearControlPanel, 'Text', 'Max shear strain \gamma_{max}:', ...
                'Position', [20 558 170 22], 'Interpreter', 'tex');
            app.GammaMaxEdit = uieditfield(app.ShearControlPanel, 'numeric', ...
                'Limits', [0.05 1.5], 'Value', 0.5, 'Position', [200 558 120 22]);
            app.PhononPlotsButton = uibutton(app.ShearControlPanel, 'push', ...
                'Text', 'Phonon Softening Plots (6 figures)', 'Position', [20 510 300 34], ...
                'FontWeight', 'bold', 'BackgroundColor', [0.47 0.67 0.19], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.runAnalysisPlots('phonon'));
            app.BondPlotsButton = uibutton(app.ShearControlPanel, 'push', ...
                'Text', 'Bond-Breaking Plots (2 figures)', 'Position', [20 466 300 34], ...
                'FontWeight', 'bold', 'BackgroundColor', [0.85 0.33 0.10], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.runAnalysisPlots('bond'));
            app.DispersionButton = uibutton(app.ShearControlPanel, 'push', ...
                'Text', '3D Dispersion Surface (Brillouin zone)', 'Position', [20 422 300 34], ...
                'FontWeight', 'bold', 'BackgroundColor', [0.20 0.45 0.80], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.openDispersionSurface());
            app.AlphaLabel = uilabel(app.ShearControlPanel, 'Text', 'Morse alpha (critical stretch 17.3%):', 'Position', [20 388 215 22]);
            app.AlphaEdit = uieditfield(app.ShearControlPanel, 'numeric', 'Limits', [1 20], 'Value', 4, 'Position', [240 388 80 22], ...
                'Tooltip', 'ASSUMED default 4. Calibrate to the real material: critical stretch = ln2/alpha (ideal tensile strain).', ...
                'ValueChangedFcn', @(s, ~) set(app.AlphaLabel, 'Text', sprintf('Morse alpha (critical stretch %.1f%%):', 100 * log(2) / s.Value)));
            app.AnalysisStatusLabel = uilabel(app.ShearControlPanel, ...
                'Text', 'Uses the lattice you built. Plots open in separate figure windows.', ...
                'Position', [20 300 310 80], 'WordWrap', 'on', 'FontColor', [0.4 0.4 0.4]);
            hdr.Tooltip = 'Sweeps shear strain on the current lattice with three bond models (harmonic, +pre-stress, Morse with bond rupture).';

            app.UIAxes3D_NoForce = uiaxes(app.ShearGrid);
            title(app.UIAxes3D_NoForce, 'Pristine Lattice (No Force)')
            app.UIAxes3D_NoForce.Layout.Row = 1;
            app.UIAxes3D_NoForce.Layout.Column = 2;

            app.UIAxes3D_WithForce = uiaxes(app.ShearGrid);
            title(app.UIAxes3D_WithForce, 'Shear-Deformed Lattice (Force Distributed)')
            app.UIAxes3D_WithForce.Layout.Row = 1;
            app.UIAxes3D_WithForce.Layout.Column = 3;

            app.UIAxesMode_NoForce = uiaxes(app.ShearGrid);
            title(app.UIAxesMode_NoForce, 'Normal Modes - Pristine')
            xlabel(app.UIAxesMode_NoForce, 'Mode Index')
            ylabel(app.UIAxesMode_NoForce, '\omega (angular frequency)')
            app.UIAxesMode_NoForce.Layout.Row = 2;
            app.UIAxesMode_NoForce.Layout.Column = 2;

            app.UIAxesMode_WithForce = uiaxes(app.ShearGrid);
            title(app.UIAxesMode_WithForce, 'Normal Modes - Sheared')
            xlabel(app.UIAxesMode_WithForce, 'Mode Index')
            ylabel(app.UIAxesMode_WithForce, '\omega (angular frequency)')
            app.UIAxesMode_WithForce.Layout.Row = 2;
            app.UIAxesMode_WithForce.Layout.Column = 3;

            app.updateShearAnalysis(struct('Value', 0));
        end

        function openDispersionSurface(app)
            if ~app.isBuilt
                app.showAlert('Please Build & Solve the lattice first.', 'Not Built');
                return;
            end
            try
                nb = size(app.AnalysisBasis, 1);
                brillouinShearSurface(app.AnalysisLV, app.AnalysisBasis, app.masses(1:nb), max(app.bonds(:, 3)), ...
                    app.BaseSpringConstant, app.GammaMaxEdit.Value, app.AlphaEdit.Value);
            catch ME
                app.showAlert(ME.message, 'Dispersion surface error');
            end
        end

        function runAnalysisPlots(app, kind)
            if ~app.isBuilt
                app.showAlert('Please Build & Solve the lattice first.', 'Not Built');
                return;
            end
            dlg = uiprogressdlg(app.ShearFig, 'Title', 'Shear analysis', ...
                'Message', 'Sweeping shear strain...', 'Indeterminate', 'on');
            cleanup = onCleanup(@() close(dlg));
            try
                gmax = app.GammaMaxEdit.Value;
                k = app.BaseSpringConstant;
                if 3 * size(app.pos, 1) > 900
                    app.showAlert(sprintf('This lattice has %d atoms; the sweeps may take several minutes. Consider a smaller N.', ...
                        size(app.pos, 1)), 'Large lattice');
                end
                c = app.AnalysisCache;
                if isempty(c) || abs(c.gammaMax - gmax) > 1e-12 || c.k ~= k || c.alphaMorse ~= app.AlphaEdit.Value
                    dlg.Message = 'Running 3 shear sweeps (affine / relaxed / Morse rupture)...';
                    c = phononShearSuite(app.pos, app.bonds, app.masses, app.AnalysisLV, app.AnalysisBasis, k, gmax, [], app.AlphaEdit.Value);
                    app.AnalysisCache = c;
                end
                if strcmp(kind, 'phonon')
                    dlg.Message = 'Drawing phonon-softening plots...';
                    plotShearSoftening(c);
                    msg = 'Phonon softening figures P1-P6 opened.';
                else
                    dlg.Message = 'Drawing bond-breaking plots...';
                    plotShearBonds(c);
                    msg = 'Bond-breaking figures Q1-Q2 opened.';
                end
                R = c.cases{c.main};
                gb = min(R.gammaFirstBreak);
                if isnan(gb), gbs = 'none'; else, gbs = sprintf('%.3f', gb); end
                if isnan(R.gammaInstab), gi = 'none'; else, gi = sprintf('%.3f', R.gammaInstab); end
                app.AnalysisStatusLabel.Text = sprintf('%s First bond breaks at gamma = %s; lowest mode goes imaginary at gamma = %s.', msg, gbs, gi);
            catch ME
                app.showAlert(ME.message, 'Analysis error');
            end
        end
        function updateShearAnalysis(app, event)
            shearVal = event.Value;

            if isempty(app.pos) || isempty(app.BaseStiffnessMatrix)
                return;
            end

            [~, elementNames] = getElementDetails(app.masses);

            k = app.BaseSpringConstant;
            [~, ~, ~, freqs_base, ~] = applyShearForce(app.pos, app.bonds, app.masses, app.BaseStiffnessMatrix, k, 0.0);
            [coords_sh, force_dist, ~, freqs_sh, ~] = applyShearForce(app.pos, app.bonds, app.masses, app.BaseStiffnessMatrix, k, shearVal);

            cla(app.UIAxes3D_NoForce);
            plotLattice(app.pos, app.bonds, app.UIAxes3D_NoForce);
            atomHover(app.UIAxes3D_NoForce, makeHoverLabels(elementNames, app.masses), @() app.pos);

            cla(app.UIAxes3D_WithForce);
            plotLatticeWithGradient(app.UIAxes3D_WithForce, coords_sh, elementNames, force_dist, app.bonds, 'Bond-force load (1 = load at gamma = 0.5)');
            atomHover(app.UIAxes3D_WithForce, makeHoverLabels(elementNames, app.masses), @() coords_sh);

            cla(app.UIAxesMode_NoForce);
            plot(app.UIAxesMode_NoForce, freqs_base, '-o', 'LineWidth', 1.5, 'Color', 'b');

            yyaxis(app.UIAxesMode_WithForce, 'left');  cla(app.UIAxesMode_WithForce);
            yyaxis(app.UIAxesMode_WithForce, 'right');  cla(app.UIAxesMode_WithForce);
            yyaxis(app.UIAxesMode_WithForce, 'left');
            app.UIAxesMode_WithForce.YColor = [0 0 0];
            hold(app.UIAxesMode_WithForce, 'on');
            hPristine = plot(app.UIAxesMode_WithForce, freqs_base, '--', 'LineWidth', 1, 'Color', [0.6 0.6 0.6]);
            hSheared = plot(app.UIAxesMode_WithForce, freqs_sh, '-o', 'LineWidth', 1.5, 'Color', 'r');
            hold(app.UIAxesMode_WithForce, 'off');
            ylabel(app.UIAxesMode_WithForce, '\omega (angular frequency)');

            yyaxis(app.UIAxesMode_WithForce, 'right');
            app.UIAxesMode_WithForce.YColor = [0.0 0.5 0.0];
            deltaOmega = freqs_sh - freqs_base;
            hold(app.UIAxesMode_WithForce, 'on');
            hDelta = plot(app.UIAxesMode_WithForce, deltaOmega, '-', 'LineWidth', 1.2, 'Color', [0.0 0.5 0.0]);
            plot(app.UIAxesMode_WithForce, [1 numel(deltaOmega)], [0 0], ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.75);
            hold(app.UIAxesMode_WithForce, 'off');
            ylabel(app.UIAxesMode_WithForce, '\Delta\omega = \omega_{sheared} - \omega_{pristine}');

            legend([hPristine hSheared hDelta], {'Pristine (dashed, left axis)', 'Sheared (left axis)', '\Delta\omega (right axis)'}, 'Location', 'northwest');

            nonzeroBase = freqs_base(freqs_base > 1e-4);
            nonzeroSh = freqs_sh(freqs_sh > 1e-4);
            if ~isempty(nonzeroBase) && ~isempty(nonzeroSh)
                lowPct = 100 * (nonzeroSh(1) - nonzeroBase(1)) / nonzeroBase(1);
                highPct = 100 * (freqs_sh(end) - freqs_base(end)) / freqs_base(end);
                title(app.UIAxesMode_WithForce, sprintf('Normal Modes - Sheared (lowest mode %+.1f%%, highest mode %+.1f%%)', ...
                    lowPct, highPct));
            end
        end

        function RelaxButtonPushed(app)
            if isempty(app.pos) || isempty(app.BaseStiffnessMatrix)
                return;
            end

            shearVal = app.ShearSlider.Value;
            k = app.BaseSpringConstant;

            try
                [coords_relaxed, coords_affine, resBefore, resAfter] = ...
                    relaxShearLattice(app.pos, app.bonds, k, app.BaseStiffnessMatrix, shearVal, 'gj');
            catch ME
                if strcmp(ME.identifier, 'gaussJordanSolve:singular')
                    app.showAlert(sprintf(['This bond network is too "floppy" to relax (not enough bonds ' ...
                        'to resist shear away from the fixed top/bottom layers - a Maxwell-rigidity ' ...
                        'issue). Rebuild with the "Auto (grows shells until rigid)" bond ' ...
                        'network and try again.\n\nDetails: %s'], ME.message), 'Cannot Relax');
                else
                    app.showAlert(ME.message, 'Relaxation error');
                end
                return;
            end

            [~, elementNames] = getElementDetails(app.masses);

            bonds_relaxed = app.bonds;
            i_idx = app.bonds(:, 1);
            j_idx = app.bonds(:, 2);
            dispField_relaxed = coords_relaxed - app.pos;
            r_old = app.bonds(:, 3) .* app.bonds(:, 4:6);
            dvec = r_old + (dispField_relaxed(j_idx, :) - dispField_relaxed(i_idx, :));
            len = sqrt(sum(dvec.^2, 2));
            valid = len > 0;
            bonds_relaxed(valid, 4:6) = dvec(valid, :) ./ len(valid);

            [K_relaxed, M_relaxed] = buildDynamicalMatrix(coords_relaxed, bonds_relaxed, k, app.masses);
            invSqrtM = diag(1 ./ sqrt(diag(M_relaxed)));
            Dmat_relaxed = invSqrtM * K_relaxed * invSqrtM;
            [~, eigVals_relaxed] = jacobiEigenSolver(Dmat_relaxed);
            freqs_relaxed = sqrt(max(sort(eigVals_relaxed), 0));

            correctionMag = sqrt(sum((coords_relaxed - coords_affine).^2, 2));
            if max(correctionMag) > 1e-9
                correctionMag = correctionMag / max(correctionMag);
            end

            cla(app.UIAxes3D_WithForce);
            plotLatticeWithGradient(app.UIAxes3D_WithForce, coords_relaxed, elementNames, correctionMag, bonds_relaxed, 'Relaxation displacement (relative)');
            atomHover(app.UIAxes3D_WithForce, makeHoverLabels(elementNames, app.masses), @() coords_relaxed);
            title(app.UIAxes3D_WithForce, 'Shear-Deformed Lattice (Relaxed, Exp.5 Gauss-Jordan)');

            yyaxis(app.UIAxesMode_WithForce, 'left');  cla(app.UIAxesMode_WithForce);
            yyaxis(app.UIAxesMode_WithForce, 'right');  cla(app.UIAxesMode_WithForce);
            yyaxis(app.UIAxesMode_WithForce, 'left');
            app.UIAxesMode_WithForce.YColor = [0 0 0];
            hold(app.UIAxesMode_WithForce, 'on');
            hPristine = plot(app.UIAxesMode_WithForce, app.omega, '--', 'LineWidth', 1, 'Color', [0.6 0.6 0.6]);
            hRelaxed = plot(app.UIAxesMode_WithForce, freqs_relaxed, '-o', 'LineWidth', 1.5, 'Color', [0.49 0.18 0.56]);
            hold(app.UIAxesMode_WithForce, 'off');
            ylabel(app.UIAxesMode_WithForce, '\omega (angular frequency)');

            yyaxis(app.UIAxesMode_WithForce, 'right');
            app.UIAxesMode_WithForce.YColor = [0.0 0.5 0.0];
            deltaOmega = freqs_relaxed - app.omega;
            hold(app.UIAxesMode_WithForce, 'on');
            hDelta = plot(app.UIAxesMode_WithForce, deltaOmega, '-', 'LineWidth', 1.2, 'Color', [0.0 0.5 0.0]);
            plot(app.UIAxesMode_WithForce, [1 numel(deltaOmega)], [0 0], ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.75);
            hold(app.UIAxesMode_WithForce, 'off');
            ylabel(app.UIAxesMode_WithForce, '\Delta\omega = \omega_{relaxed} - \omega_{pristine}');

            legend([hPristine hRelaxed hDelta], {'Pristine (dashed, left axis)', 'Sheared, relaxed (left axis)', '\Delta\omega (right axis)'}, 'Location', 'northwest');

            nonzeroBase = app.omega(app.omega > 1e-4);
            nonzeroRelaxed = freqs_relaxed(freqs_relaxed > 1e-4);
            if ~isempty(nonzeroBase) && ~isempty(nonzeroRelaxed)
                lowPct = 100 * (nonzeroRelaxed(1) - nonzeroBase(1)) / nonzeroBase(1);
                highPct = 100 * (freqs_relaxed(end) - app.omega(end)) / app.omega(end);
                title(app.UIAxesMode_WithForce, sprintf('Normal Modes - Sheared (Relaxed) (lowest mode %+.1f%%, highest mode %+.1f%%)', ...
                    lowPct, highPct));
            else
                title(app.UIAxesMode_WithForce, 'Normal Modes - Sheared (Relaxed)');
            end

            app.RelaxStatusLabel.Text = sprintf(['Max residual force per atom: %.3e before relax -> ' ...
                '%.3e after (Gauss-Jordan). Move the slider to go back to the live affine view.'], ...
                resBefore, resAfter);
            app.RelaxStatusLabel.FontColor = [0.1 0.5 0.1];
        end

        function switchToView(app, which)
            app.ensurePlotFigure();
            is3D = strcmp(which, '3D');

            if is3D
                app.Panel3D.Visible = 'on';
                app.PanelSpec.Visible = 'off';
            else
                app.Panel3D.Visible = 'off';
                app.PanelSpec.Visible = 'on';
            end

            activeColor = [0.20 0.45 0.80];
            activeFont = 'w';
            inactiveColor = [0.94 0.94 0.94];
            inactiveFont = 'k';

            if is3D
                app.View3DButton.BackgroundColor = activeColor;
                app.View3DButton.ForegroundColor = activeFont;
                app.ViewSpecButton.BackgroundColor = inactiveColor;
                app.ViewSpecButton.ForegroundColor = inactiveFont;
            else
                app.View3DButton.BackgroundColor = inactiveColor;
                app.View3DButton.ForegroundColor = inactiveFont;
                app.ViewSpecButton.BackgroundColor = activeColor;
                app.ViewSpecButton.ForegroundColor = activeFont;
            end

            figure(app.PlotFig);
        end

        function LoadFileButtonPushed(app)
            [file, path] = uigetfile( ...
                {'*.txt;*.dat;*.cif', 'Crystal files (*.txt, *.dat, *.cif)'; ...
                 '*.cif', 'CIF files (*.cif)'; ...
                 '*.txt;*.dat;*.cif', 'Custom crystal unit files (*.txt, *.dat,*.cif)'}, ...
                'Select a crystal unit file');
            if isequal(file, 0)
                return;
            end
            fullPath = fullfile(path, file);
            [~, ~, ext] = fileparts(file);
            try
                if strcmpi(ext, '.cif')
                    [app.latticeVectors, app.basisFrac, app.basisMasses] = parseCIFFile(fullPath);
                else
                    [app.latticeVectors, app.basisFrac, app.basisMasses] = parseCrystalFile(fullPath);
                end
                app.useCustomFile = true;
                app.StatusLabel.Text = sprintf('Loaded: %s (%d atoms/cell)', ...
                    file, size(app.basisFrac, 1));
                app.StatusLabel.FontColor = [0.1 0.5 0.1];
            catch ME
                app.showAlert(ME.message, 'Could not load file');
            end
        end

        function BuildButtonPushed(app)
            app.stopTimerIfRunning();
            try
                N = round(app.NEditField.Value);
                k = app.kEditField.Value;

                if app.useCustomFile
                    latVec = app.latticeVectors;
                    basis = app.basisFrac;
                    bMasses = app.basisMasses;
                else
                    [latVec, basis, bMasses] = getBuiltinCrystalDef('bcc');
                end

                app.AnalysisLV = latVec;
                app.AnalysisBasis = basis;
                app.AnalysisCache = [];

                shellsUsed = [];
                if strcmp(app.cutoffDropDown.Value, 'Auto (grows shells until rigid)')
                    app.setStatus('Building lattice (auto-growing shells until rigid)...', [0.2 0.2 0.6]);
                    [app.pos, app.bonds, app.masses, shellsUsed] = autoRigidShells(latVec, basis, bMasses, N, k);
                else
                    if strcmp(app.cutoffDropDown.Value, 'Nearest neighbors only')
                        numShells = 1;
                    else
                        numShells = 2;
                    end
                    app.setStatus('Building lattice...', [0.2 0.2 0.6]);
                    [app.pos, app.bonds, app.masses] = generateLatticeGeneral(latVec, basis, bMasses, N, numShells);
                end

                app.setStatus('Assembling dynamical matrix...', [0.2 0.2 0.6]);
                [K, ~, Dmat] = buildDynamicalMatrix(app.pos, app.bonds, k, app.masses);
                app.BaseStiffnessMatrix = K;
                app.BaseSpringConstant = k;

                app.setStatus('Solving eigenproblem...', [0.2 0.2 0.6]);
                [V, eigVals] = jacobiEigenSolver(Dmat, 1e-12, 200);
                [eigVals, order] = sort(eigVals);
                app.V = V(:, order);
                app.omega = sqrt(max(eigVals, 0));

                app.isBuilt = true;
                app.populateModeDropdown();
                app.plotStaticLattice();
                app.plotFrequencySpectrum();
                app.ModeDropDownValueChanged();

                Natoms = size(app.pos, 1);
                numZero = sum(eigVals < 1e-6);
                if isempty(shellsUsed)
                    app.setStatus(sprintf('Done: %d atoms, %d bonds, %d modes (%d rigid-body zero modes).', ...
                        Natoms, size(app.bonds, 1), length(app.omega), numZero), [0.1 0.5 0.1]);
                else
                    app.setStatus(sprintf('Done: %d atoms, %d bonds, %d modes (%d rigid-body zero modes, auto used %d shell(s)).', ...
                        Natoms, size(app.bonds, 1), length(app.omega), numZero, shellsUsed), [0.1 0.5 0.1]);
                end

                app.ModeDropDown.Enable = 'on';
                app.PlayButton.Enable = 'on';
                app.PauseButton.Enable = 'off';
                app.StopButton.Enable = 'off';
                app.LaunchShearButton.Enable = 'on';
            catch ME
                app.setStatus('Build failed - see error dialog.', [0.7 0.1 0.1]);
                app.showAlert(ME.message, 'Build error');
            end
        end

        function populateModeDropdown(app)
            nModes = length(app.omega);
            items = cell(1, nModes);
            for idx = 1:nModes
                tag = '';
                if app.omega(idx) < 1e-6
                    tag = ' (rigid body)';
                end
                items{idx} = sprintf('Mode %d - omega = %.4f%s', idx, app.omega(idx), tag);
            end

            app.ModeDropDown.Items = items;

            firstVibrating = find(app.omega > 1e-6, 1);
            if isempty(firstVibrating)
                firstVibrating = 1;
            end
            app.ModeDropDown.Value = items{firstVibrating};
        end

        function modeIdx = getSelectedModeIndex(app)
            val = app.ModeDropDown.Value;
            if isempty(val)
                modeIdx = 1;
                return;
            end
            if isnumeric(val)
                modeIdx = val;
            else
                modeIdx = find(strcmp(app.ModeDropDown.Items, val), 1);
                if isempty(modeIdx)
                    modeIdx = 1;
                end
            end
        end

        function plotFrequencySpectrum(app)
            app.ensurePlotFigure();
            cla(app.SpecAx);
            hold(app.SpecAx, 'on');

            nModes = length(app.omega);
            idxAll = 1:nModes;
            plot(app.SpecAx, idxAll, app.omega, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor', [0.2 0.45 0.8]);
            hold(app.SpecAx, 'off');
            app.SpecHighlightHandle = [];
        end

        function plotStaticLattice(app)
            app.ensurePlotFigure();
            ax = app.PlotAx;
            cla(ax);
            hold(ax, 'on');

            [elementColors, elementNames] = getElementDetails(app.masses);

            nBonds = size(app.bonds, 1);
            app.bondHandles = gobjects(nBonds, 1);
            for b = 1:nBonds
                i = app.bonds(b, 1);
                j = app.bonds(b, 2);
                app.bondHandles(b) = plot3(ax, ...
                    [app.pos(i,1) app.pos(j,1)], ...
                    [app.pos(i,2) app.pos(j,2)], ...
                    [app.pos(i,3) app.pos(j,3)], ...
                    '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
            end

            app.atomHandle = scatter3(ax, app.pos(:,1), app.pos(:,2), app.pos(:,3), ...
                160, elementColors, 'filled', 'MarkerEdgeColor', 'k');

            hoverTxt = makeHoverLabels(elementNames, app.masses);
            ah = app.atomHandle;
            atomHover(ax, hoverTxt, @() [ah.XData(:), ah.YData(:), ah.ZData(:)]);

            axis(ax, 'equal');
            grid(ax, 'on');
            view(ax, 3);
            title(ax, sprintf('%d atoms, %d bonds (hover an atom for details)', size(app.pos, 1), nBonds));
            hold(ax, 'off');
        end

        function highlightModeOnSpectrum(app, modeIdx)
            app.ensurePlotFigure();
            if isempty(app.SpecHighlightHandle) || ~isvalid(app.SpecHighlightHandle)
                hold(app.SpecAx, 'on');
                app.SpecHighlightHandle = plot(app.SpecAx, modeIdx, app.omega(modeIdx), ...
                    'o', 'MarkerSize', 12, 'MarkerFaceColor', [0.85 0.2 0.2], ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
                hold(app.SpecAx, 'off');
            else
                set(app.SpecHighlightHandle, 'XData', modeIdx, 'YData', app.omega(modeIdx));
            end
        end

        function ModeDropDownValueChanged(app)
            if ~app.isBuilt
                return;
            end
            modeIdx = app.getSelectedModeIndex();
            Natoms = size(app.pos, 1);
            app.modeShape = real(reshape(app.V(:, modeIdx), 3, Natoms))';
            app.currentOmega = app.omega(modeIdx);
            app.animTime = 0;
            app.highlightModeOnSpectrum(modeIdx);
        end

        function PlayButtonPushed(app)
            if ~app.isBuilt
                return;
            end

            if isempty(app.modeShape)
                app.ModeDropDownValueChanged();
            end

            if isempty(app.animTimer) || ~isvalid(app.animTimer)
                app.animTimer = timer( ...
                    'ExecutionMode', 'fixedRate', ...
                    'BusyMode', 'drop', ...
                    'Period', app.speedToPeriod(), ...
                    'TimerFcn', @(~, ~) app.animationTick());
            end
            start(app.animTimer);

            app.PlayButton.Enable = 'off';
            app.PauseButton.Enable = 'on';
            app.StopButton.Enable = 'on';
        end

        function PauseButtonPushed(app)
            if ~isempty(app.animTimer) && isvalid(app.animTimer) && strcmp(app.animTimer.Running, 'on')
                stop(app.animTimer);
            end
            app.PlayButton.Enable = 'on';
            app.PauseButton.Enable = 'off';
        end

        function StopButtonPushed(app)
            app.stopTimerIfRunning();
            app.animTime = 0;
            if app.isBuilt
                app.plotStaticLattice();
            end
            app.PlayButton.Enable = 'on';
            app.PauseButton.Enable = 'off';
            app.StopButton.Enable = 'off';
        end

        function SpeedSliderValueChanged(app)
        end

        function period = speedToPeriod(app)
            period = 1 / app.AnimFrameRateHz;
        end

        function animationTick(app)
            if isempty(app.modeShape) || isempty(app.atomHandle) || ~isvalid(app.atomHandle)
                return;
            end

            displayOmega = 2 * pi * app.SpeedSlider.Value;
            app.animTime = app.animTime + app.speedToPeriod();
            amp = app.AmplitudeSlider.Value;
            displaced = app.pos + amp * app.modeShape * cos(displayOmega * app.animTime);

            set(app.atomHandle, 'XData', displaced(:,1), 'YData', displaced(:,2), 'ZData', displaced(:,3));

            for b = 1:numel(app.bondHandles)
                if isvalid(app.bondHandles(b))
                    i = app.bonds(b, 1);
                    j = app.bonds(b, 2);
                    set(app.bondHandles(b), ...
                        'XData', [displaced(i,1) displaced(j,1)], ...
                        'YData', [displaced(i,2) displaced(j,2)], ...
                        'ZData', [displaced(i,3) displaced(j,3)]);
                end
            end
            drawnow limitrate;
        end

        function stopTimerIfRunning(app)
            if ~isempty(app.animTimer) && isvalid(app.animTimer)
                stop(app.animTimer);
                delete(app.animTimer);
            end
            app.animTimer = [];
        end

        function setStatus(app, msg, color)
            if nargin < 3
                color = [0.4 0.4 0.4];
            end
            app.StatusLabel.Text = msg;
            app.StatusLabel.FontColor = color;
            drawnow;
        end

        function showAlert(app, msg, titleStr)
            if ~isempty(app.UIFigure) && isvalid(app.UIFigure)
                uialert(app.UIFigure, msg, titleStr);
            else
                errordlg(msg, titleStr);
            end
        end

        function onClose(app)
            app.delete();
        end
    end
end

function labels = makeHoverLabels(elementNames, masses)
    n = numel(masses);
    labels = cell(n, 1);
    for a = 1:n
        nm = elementNames{a};
        tok = regexp(nm, '^(.*)\s\((\w+)\)$', 'tokens', 'once');
        if isempty(tok)
            labels{a} = sprintf('%s\natom #%d', nm, a);
        else
            labels{a} = sprintf('%s  -  %s\nmass %.3f amu   atom #%d', tok{2}, tok{1}, masses(a), a);
        end
    end
end
