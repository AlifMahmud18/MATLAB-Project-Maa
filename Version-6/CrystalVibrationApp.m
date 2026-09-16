classdef CrystalVibrationApp < handle
    %CRYSTALVIBRATIONAPP Interactive crystal vibration viewer with element colors, hover tooltips, and shear analysis.

    properties (Access = private)
        % --- Main UI components ---
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

        % --- Combined plot window (3D view + frequency spectrum) ---
        PlotFig
        Panel3D
        PanelSpec
        PlotAx
        SpecAx
        View3DButton
        ViewSpecButton
        SpecHighlightHandle

        % --- Shear Analysis UI components ---
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

        % --- Crystal definition currently loaded ---
        latticeVectors
        basisFrac
        basisMasses
        useCustomFile = false

        % --- Solved model ---
        pos
        bonds
        masses         
        V              
        omega          
        isBuilt = false
        BaseStiffnessMatrix % Stored for shear analysis
        BaseSpringConstant  % k used to build BaseStiffnessMatrix, kept in sync for shear analysis

        % --- Animation state ---
        modeShape
        currentOmega = 0
        animTime = 0
        animTimer
        atomHandle
        bondHandles
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

        %% ---------------- UI construction ----------------
        function createComponents(app)
            app.UIFigure = uifigure('Name', 'Crystal Vibration Explorer - Controls', ...
                'Position', [100 100 320 650], 'Color', 'w');
            app.UIFigure.CloseRequestFcn = @(~, ~) app.onClose();

            uilabel(app.UIFigure, 'Text', 'Crystal Vibration Explorer', ...
                'FontSize', 16, 'FontWeight', 'bold', 'Position', [20 610 280 25]);

            % --- Load file section ---
            app.LoadFileButton = uibutton(app.UIFigure, 'push', ...
                'Text', 'Load Crystal File...', 'Position', [20 575 260 30], ...
                'ButtonPushedFcn', @(~, ~) app.LoadFileButtonPushed());

            app.StatusLabel = uilabel(app.UIFigure, 'Text', 'No file loaded (using default bcc structure)', ...
                'Position', [20 550 280 20], 'FontColor', [0.4 0.4 0.4]);

            % --- Physics parameters ---
            uilabel(app.UIFigure, 'Text', 'Supercell size N (NxNxN cells):', 'Position', [20 515 260 20]);
            app.NEditField = uieditfield(app.UIFigure, 'numeric', ...
                'Value', 2, 'Limits', [1 4], 'Position', [20 490 260 25]);

            uilabel(app.UIFigure, 'Text', 'Spring constant k:', 'Position', [20 460 260 20]);
            app.kEditField = uieditfield(app.UIFigure, 'numeric', ...
                'Value', 1.0, 'Limits', [0.001 Inf], 'Position', [20 435 260 25]);

            uilabel(app.UIFigure, 'Text', 'Bond network:', 'Position', [20 400 260 20]);
            app.cutoffDropDown = uidropdown(app.UIFigure, ...
                'Items', {'Nearest neighbors only', '1st + 2nd neighbors (recommended, rigid)'}, ...
                'Value', '1st + 2nd neighbors (recommended, rigid)', ...
                'Position', [20 375 260 25]);

            app.BuildButton = uibutton(app.UIFigure, 'push', ...
                'Text', 'Build && Solve', 'Position', [20 330 260 35], ...
                'FontWeight', 'bold', 'BackgroundColor', [0.20 0.45 0.80], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.BuildButtonPushed());

            % --- Mode selection ---
            uilabel(app.UIFigure, 'Text', 'Vibrational mode:', 'Position', [20 295 260 20]);
            app.ModeDropDown = uidropdown(app.UIFigure, ...
                'Items', {'(build the model first)'}, 'Enable', 'off', ...
                'Position', [20 270 260 25], ...
                'ValueChangedFcn', @(~, ~) app.ModeDropDownValueChanged());

            % --- Animation controls ---
            uilabel(app.UIFigure, 'Text', 'Amplitude:', 'Position', [20 235 260 20]);
            app.AmplitudeSlider = uislider(app.UIFigure, ...
                'Limits', [0 1], 'Value', 0.25, 'Position', [25 220 250 3]);

            uilabel(app.UIFigure, 'Text', 'Speed (fps):', 'Position', [20 175 260 20]);
            app.SpeedSlider = uislider(app.UIFigure, ...
                'Limits', [5 60], 'Value', 30, 'Position', [25 160 250 3], ...
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

            % --- Shear Analysis Launcher ---
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

        %% ---------------- Shear Window Construction ----------------
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
            app.ShearLabel.Position = [20 750 150 22];
            app.ShearLabel.Text = 'Shear Force Magnitude:';

            app.ShearSlider = uislider(app.ShearControlPanel);
            app.ShearSlider.Limits = [0 2];
            app.ShearSlider.Position = [20 720 300 3];
            app.ShearSlider.ValueChangedFcn = @(src, event) app.updateShearAnalysis(event);

            % --- Relaxation (Exp. 5: Gauss-Jordan / Gauss-Seidel) ---
            % The slider above only ever shows the purely affine (Cauchy-
            % Born) guess, recomputed live on every drag. Relaxation is a
            % one-shot linear solve (relaxShearLattice.m) instead - too
            % slow to run on every slider tick for a larger supercell, so
            % it is a separate, on-demand button rather than being wired
            % into ValueChangedFcn.
            app.RelaxButton = uibutton(app.ShearControlPanel, 'push', ...
                'Text', 'Relax (Gauss-Jordan, Exp.5)', 'Position', [20 670 300 32], ...
                'FontWeight', 'bold', 'BackgroundColor', [0.20 0.45 0.80], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.RelaxButtonPushed());

            app.RelaxStatusLabel = uilabel(app.ShearControlPanel);
            app.RelaxStatusLabel.Position = [20 630 310 34];
            app.RelaxStatusLabel.Text = 'Relax to solve for true equilibrium of interior atoms.';
            app.RelaxStatusLabel.FontColor = [0.4 0.4 0.4];
            app.RelaxStatusLabel.WordWrap = 'on';

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

        function updateShearAnalysis(app, event)
            shearVal = event.Value;
            
            if isempty(app.pos) || isempty(app.BaseStiffnessMatrix)
                return;
            end
            
            % Assume element types array based on masses (used for coloring, not physics)
            [~, elementNames] = getElementDetails(app.masses);
            
            % applyShearForce indexes K_matrix in 3x3 per-atom blocks, so it
            % needs the assembled stiffness/dynamical matrix, not a scalar.
            % It also needs the bond list and the same scalar spring
            % constant used to build that matrix, so it can rebuild the
            % dynamical matrix for the sheared geometry.
            k = app.BaseSpringConstant;
            [~, ~, ~, freqs_base, ~] = applyShearForce(app.pos, app.bonds, app.masses, app.BaseStiffnessMatrix, k, 0.0);
            [coords_sh, force_dist, ~, freqs_sh, ~] = applyShearForce(app.pos, app.bonds, app.masses, app.BaseStiffnessMatrix, k, shearVal);
            
            cla(app.UIAxes3D_NoForce);
            plotLattice(app.pos, app.bonds, app.UIAxes3D_NoForce);
            
            cla(app.UIAxes3D_WithForce);
            plotLatticeWithGradient(app.UIAxes3D_WithForce, coords_sh, elementNames, force_dist);
            
            cla(app.UIAxesMode_NoForce);
            plot(app.UIAxesMode_NoForce, freqs_base, '-o', 'LineWidth', 1.5, 'Color', 'b');
            
            cla(app.UIAxesMode_WithForce);
            plot(app.UIAxesMode_WithForce, freqs_sh, '-o', 'LineWidth', 1.5, 'Color', 'r');
        end

        %% ---------------- Relaxation (Exp. 5 linear solve) ----------------
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
                        'issue). Rebuild with the "1st + 2nd neighbors (recommended, rigid)" bond ' ...
                        'network and try again.\n\nDetails: %s'], ME.message), 'Cannot Relax');
                else
                    app.showAlert(ME.message, 'Relaxation error');
                end
                return;
            end

            [~, elementNames] = getElementDetails(app.masses);

            % Rebuild bond unit vectors for the RELAXED geometry, then
            % rebuild K/Dmat and re-solve for the relaxed vibrational
            % spectrum.
            %
            % This must NOT be done via
            % "coords_relaxed(j,:) - coords_relaxed(i,:)" (the same
            % wraparound pitfall documented and fixed in
            % computeShearForces.m / applyShearForce.m): bonds can
            % connect an atom to a periodic image of its neighbor, so
            % the correct new bond vector is the ORIGINAL (minimum-
            % image) bond vector plus the per-atom displacement
            % DIFFERENCE (a well-posed quantity here, since every atom's
            % displacement - affine for boundary atoms, solved for free
            % atoms - is already self-consistent with the same K used to
            % find it).
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

            % Color by the SIZE of the relaxation correction itself
            % (|coords_relaxed - coords_affine|) rather than residual
            % force: free atoms are solved to have ~zero net force by
            % construction, so a force-colored plot would just look
            % uniformly "settled" - the correction magnitude instead
            % shows WHERE the affine guess was wrong and by how much.
            correctionMag = sqrt(sum((coords_relaxed - coords_affine).^2, 2));
            if max(correctionMag) > 0
                correctionMag = correctionMag / max(correctionMag);
            end

            cla(app.UIAxes3D_WithForce);
            plotLatticeWithGradient(app.UIAxes3D_WithForce, coords_relaxed, elementNames, correctionMag);
            title(app.UIAxes3D_WithForce, 'Shear-Deformed Lattice (Relaxed, Exp.5 Gauss-Jordan)');

            cla(app.UIAxesMode_WithForce);
            plot(app.UIAxesMode_WithForce, freqs_relaxed, '-o', 'LineWidth', 1.5, 'Color', [0.49 0.18 0.56]);
            title(app.UIAxesMode_WithForce, 'Normal Modes - Sheared (Relaxed)');

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

        %% ---------------- Callbacks: loading a crystal ----------------
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

        %% ---------------- Callback: build & solve ----------------
        function BuildButtonPushed(app)
            app.stopTimerIfRunning();
            try
                N = round(app.NEditField.Value);
                k = app.kEditField.Value;

                if strcmp(app.cutoffDropDown.Value, 'Nearest neighbors only')
                    numShells = 1;
                else
                    numShells = 2;
                end

                if app.useCustomFile
                    latVec = app.latticeVectors;
                    basis = app.basisFrac;
                    bMasses = app.basisMasses;
                else
                    [latVec, basis, bMasses] = getBuiltinCrystalDef('bcc');
                end

                app.setStatus('Building lattice...', [0.2 0.2 0.6]);
                [app.pos, app.bonds, app.masses] = generateLatticeGeneral(latVec, basis, bMasses, N, numShells);

                app.setStatus('Assembling dynamical matrix...', [0.2 0.2 0.6]);
                [~, ~, Dmat] = buildDynamicalMatrix(app.pos, app.bonds, k, app.masses);
                app.BaseStiffnessMatrix = Dmat; % Saved for shear matrix analysis
                app.BaseSpringConstant = k;     % kept alongside Dmat so shear analysis stays consistent even if kEditField changes later

                app.setStatus('Solving eigenproblem...', [0.2 0.2 0.6]);
                [V, eigVals] = jacobiEigenSolver(Dmat, 1e-12, 200);
                [eigVals, order] = sort(eigVals);
                app.V = V(:, order);
                app.omega = sqrt(max(eigVals, 0));

                app.isBuilt = true;
                app.populateModeDropdown();
                app.plotStaticLattice();
                app.plotFrequencySpectrum();
                app.ModeDropDownValueChanged(); % prime modeShape/currentOmega + spectrum highlight

                Natoms = size(app.pos, 1);
                numZero = sum(eigVals < 1e-6);
                app.setStatus(sprintf('Done: %d atoms, %d bonds, %d modes (%d rigid-body zero modes).', ...
                    Natoms, size(app.bonds, 1), length(app.omega), numZero), [0.1 0.5 0.1]);

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
            app.SpecHighlightHandle = []; % axis was cleared, so any old marker handle is stale
        end

        %% ---------------- Static 3D lattice plot (with element colors + hover tooltips) ----------------
        function plotStaticLattice(app)
            app.ensurePlotFigure();
            ax = app.PlotAx;
            cla(ax);
            hold(ax, 'on');

            [elementColors, elementNames] = getElementDetails(app.masses);

            % Bonds first so atom markers render on top of them
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

            % Hover tooltips: element, mass, and atom index
            app.atomHandle.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Element', elementNames);
            app.atomHandle.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Mass', app.masses);
            app.atomHandle.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Atom #', (1:numel(app.masses))');
            datacursormode(app.PlotFig, 'on');

            axis(ax, 'equal');
            grid(ax, 'on');
            view(ax, 3);
            title(ax, sprintf('%d atoms, %d bonds (hover an atom for details)', size(app.pos, 1), nBonds));
            hold(ax, 'off');
        end

        %% ---------------- Mode spectrum highlight ----------------
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

        %% ---------------- Mode selection ----------------
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

        %% ---------------- Animation controls ----------------
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
                app.plotStaticLattice(); % snap back to rest positions
            end
            app.PlayButton.Enable = 'on';
            app.PauseButton.Enable = 'off';
            app.StopButton.Enable = 'off';
        end

        function SpeedSliderValueChanged(app)
            if ~isempty(app.animTimer) && isvalid(app.animTimer)
                wasRunning = strcmp(app.animTimer.Running, 'on');
                stop(app.animTimer);
                app.animTimer.Period = app.speedToPeriod();
                if wasRunning
                    start(app.animTimer);
                end
            end
        end

        function period = speedToPeriod(app)
            fps = max(app.SpeedSlider.Value, 1);
            period = round((1 / fps) * 1000) / 1000; % timer periods must be multiples of 1 ms
            period = max(period, 0.001);
        end

        function animationTick(app)
            if isempty(app.modeShape) || isempty(app.atomHandle) || ~isvalid(app.atomHandle)
                return;
            end

            app.animTime = app.animTime + app.speedToPeriod();
            amp = app.AmplitudeSlider.Value;
            displaced = app.pos + amp * app.modeShape * cos(app.currentOmega * app.animTime);

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

        %% ---------------- Status / alerts / window close ----------------
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

%% ==================== Local helper functions ====================
% getElementDetails is a local function (not a class method), which is
% why it's called as a plain function (getElementDetails(...)) rather
% than app.getElementDetails(...).

function [elementColors, elementNames] = getElementDetails(masses)
%GETELEMENTDETAILS Guess an element identity/display color for each atom
%   from its mass (amu), using nearest match against a small reference
%   table. Used to color atoms in the 3D views and populate hover tips.

    table = { ...
        'H',   1.008,  [1.00 1.00 1.00]; ...
        'C',  12.011,  [0.30 0.30 0.30]; ...
        'N',  14.007,  [0.20 0.20 0.90]; ...
        'O',  15.999,  [0.90 0.10 0.10]; ...
        'Na', 22.990,  [0.60 0.20 0.80]; ...
        'Mg', 24.305,  [0.20 0.70 0.20]; ...
        'Al', 26.982,  [0.60 0.60 0.65]; ...
        'Si', 28.085,  [0.80 0.70 0.40]; ...
        'Cl', 35.450,  [0.10 0.80 0.10]; ...
        'K',  39.098,  [0.55 0.15 0.70]; ...
        'Ca', 40.078,  [0.35 0.55 0.35]; ...
        'Fe', 55.845,  [0.72 0.45 0.20]; ...
        'Cu', 63.546,  [0.72 0.30 0.10]; ...
        'Zn', 65.380,  [0.40 0.45 0.50]; ...
        'Ag', 107.868, [0.75 0.75 0.75]; ...
        'Au', 196.967, [0.85 0.65 0.13]; ...
        'Pb', 207.200, [0.30 0.30 0.35] ...
        };

    refMass = cell2mat(table(:, 2));
    n = numel(masses);
    elementNames = cell(n, 1);
    elementColors = zeros(n, 3);

    for a = 1:n
        [~, k] = min(abs(refMass - masses(a)));
        elementNames{a} = table{k, 1};
        elementColors(a, :) = table{k, 3};
    end
end

% applyShearForce, plotLattice, and plotLatticeWithGradient are provided
% as separate files (applyShearForce.m, plotLattice.m,
% plotLatticeWithGradient.m) on the MATLAB path. They are intentionally
% NOT redefined here as local functions: a local function of the same
% name in this class file would shadow the real path versions for every
% call made from within this class, silently substituting the wrong
% implementation. Make sure those three files sit next to this one (or
% anywhere else on the path) alongside buildDynamicalMatrix.m,
% jacobiEigenSolver.m, parseCIFFile.m, parseCrystalFile.m, and
% getBuiltinCrystalDef.m, none of which are defined in this file either.