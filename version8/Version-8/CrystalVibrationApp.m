classdef CrystalVibrationApp < handle
    %CRYSTALVIBRATIONAPP Single-window GUI: load a CIF/txt crystal, build the spring network, solve the modes and animate them.
    %   The whole interface lives in one uifigure - a control column on the left and a tab group on the right.
    %   Every component is placed by uigridlayout, so nothing overlaps at any window size.
    %   Build & Solve: autoRigidShells/generateLatticeGeneral -> buildDynamicalMatrix -> jacobiEigenSolver.
    %   "Lattice & Modes" tab: the 3D lattice (hover an atom for its element, via atomHover + getElementDetails)
    %   and the frequency spectrum, both visible at the same time.
    %   "Shear Analysis" tab: shear slider (applyShearForce), Relax button (relaxShearLattice), the pristine and
    %   sheared lattices with their normal-mode spectra, and the V7 buttons that run phononShearSuite.
    %   Those buttons open no windows either: "Phonon Softening" (P1-P6), "Bond Breaking" (Q1-Q2) and
    %   "Dispersion Surface" are tabs of this same figure, filled by plotShearSoftening / plotShearBonds /
    %   brillouinShearSurface drawing into the tab they are handed.

    properties (Access = private)
        UIFigure
        TabGroup
        LatticeTab
        ShearTab
        PhononTab
        BondTab
        DispersionTab
        StatusLabel

        % Left control column
        LoadFileButton
        NEditField
        kEditField
        cutoffDropDown
        BuildButton
        ModeDropDown
        AmplitudeSlider
        SpeedSlider
        PlayPauseButton
        StopButton

        % "Lattice & Modes" tab
        PlotAx
        SpecAx
        SpecHighlightHandle

        % "Shear Analysis" tab
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
        shearViewDirty = true
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
            if ~isempty(app.UIFigure) && isvalid(app.UIFigure)
                delete(app.UIFigure);
            end
        end
    end

    methods (Access = private)

        function createComponents(app)
            app.UIFigure = uifigure('Name', 'Crystal Vibration Explorer', ...
                'Position', [60 60 1500 880], 'Color', 'w');
            app.UIFigure.CloseRequestFcn = @(~, ~) app.onClose();

            mainGrid = uigridlayout(app.UIFigure, [2 2]);
            mainGrid.ColumnWidth = {340, '1x'};
            mainGrid.RowHeight = {'1x', 24};
            mainGrid.RowSpacing = 6;
            mainGrid.ColumnSpacing = 8;

            controlPanel = uipanel(mainGrid, 'Title', 'Model & Animation', ...
                'BackgroundColor', 'w');
            controlPanel.Layout.Row = 1;
            controlPanel.Layout.Column = 1;
            app.buildControlColumn(controlPanel);

            app.TabGroup = uitabgroup(mainGrid);
            app.TabGroup.Layout.Row = 1;
            app.TabGroup.Layout.Column = 2;
            app.TabGroup.SelectionChangedFcn = @(~, event) app.TabSelectionChanged(event);

            app.LatticeTab = uitab(app.TabGroup, 'Title', 'Lattice & Modes');
            app.buildLatticeTab(app.LatticeTab);

            app.ShearTab = uitab(app.TabGroup, 'Title', 'Shear Analysis');
            app.buildShearTab(app.ShearTab);

            app.PhononTab = uitab(app.TabGroup, 'Title', 'Phonon Softening', 'BackgroundColor', 'w');
            app.BondTab = uitab(app.TabGroup, 'Title', 'Bond Breaking', 'BackgroundColor', 'w');
            app.DispersionTab = uitab(app.TabGroup, 'Title', 'Dispersion Surface', 'BackgroundColor', 'w');
            app.resetAnalysisTabs();

            app.StatusLabel = uilabel(mainGrid, ...
                'Text', 'No file loaded (using default bcc structure).', ...
                'FontColor', [0.4 0.4 0.4]);
            app.StatusLabel.Layout.Row = 2;
            app.StatusLabel.Layout.Column = [1 2];

            app.setShearControlsEnabled(false);
        end

        function buildControlColumn(app, parent)
            g = uigridlayout(parent, [18 1]);
            g.ColumnWidth = {'1x'};
            g.RowHeight = {24, 32, 22, 24, 22, 24, 22, 24, 38, 26, 22, 24, 22, 50, 22, 50, 34, '1x'};
            g.RowSpacing = 6;
            g.Padding = [10 10 10 10];
            g.Scrollable = 'on';

            lbl = uilabel(g, 'Text', 'Crystal & bond network', 'FontWeight', 'bold');
            lbl.Layout.Row = 1;

            app.LoadFileButton = uibutton(g, 'push', 'Text', 'Load Crystal File...', ...
                'ButtonPushedFcn', @(~, ~) app.LoadFileButtonPushed());
            app.LoadFileButton.Layout.Row = 2;

            lbl = uilabel(g, 'Text', 'Supercell size (fixed): 1 x 1 x 1 cell');
            lbl.Layout.Row = 3;
            app.NEditField = uieditfield(g, 'numeric', 'Value', 1, 'Limits', [1 1], 'Editable', 'off');
            app.NEditField.Layout.Row = 4;

            lbl = uilabel(g, 'Text', 'Spring constant k:');
            lbl.Layout.Row = 5;
            % Spring constant k is a reduced-unit default (1.0): use a force constant fitted to the real material.
            app.kEditField = uieditfield(g, 'numeric', 'Value', 1.0, 'Limits', [0.001 Inf]);
            app.kEditField.Layout.Row = 6;

            lbl = uilabel(g, 'Text', 'Bond network:');
            lbl.Layout.Row = 7;
            app.cutoffDropDown = uidropdown(g, ...
                'Items', {'Nearest neighbors only', '1st + 2nd neighbors (fixed)', 'Auto (grows shells until rigid)'}, ...
                'Value', 'Auto (grows shells until rigid)');
            app.cutoffDropDown.Layout.Row = 8;

            app.BuildButton = uibutton(g, 'push', 'Text', 'Build && Solve', ...
                'FontWeight', 'bold', 'BackgroundColor', [0.20 0.45 0.80], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.BuildButtonPushed());
            app.BuildButton.Layout.Row = 9;

            lbl = uilabel(g, 'Text', 'Mode animation', 'FontWeight', 'bold');
            lbl.Layout.Row = 10;

            lbl = uilabel(g, 'Text', 'Vibrational mode:');
            lbl.Layout.Row = 11;
            app.ModeDropDown = uidropdown(g, 'Items', {'(build the model first)'}, 'Enable', 'off', ...
                'ValueChangedFcn', @(~, ~) app.ModeDropDownValueChanged());
            app.ModeDropDown.Layout.Row = 12;

            lbl = uilabel(g, 'Text', 'Amplitude:');
            lbl.Layout.Row = 13;
            app.AmplitudeSlider = uislider(g, 'Limits', [0 1], 'Value', 0.25);
            app.AmplitudeSlider.Layout.Row = 14;

            lbl = uilabel(g, 'Text', 'Playback speed (cycles/sec):');
            lbl.Layout.Row = 15;
            app.SpeedSlider = uislider(g, 'Limits', [0.1 10], 'Value', 2);
            app.SpeedSlider.Layout.Row = 16;

            transport = uigridlayout(g, [1 2]);
            transport.Layout.Row = 17;
            transport.ColumnWidth = {'1x', '1x'};
            transport.ColumnSpacing = 8;
            transport.Padding = [0 0 0 0];

            app.PlayPauseButton = uibutton(transport, 'push', 'Text', 'Play', 'Enable', 'off', ...
                'BackgroundColor', [0.30 0.70 0.35], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.PlayPauseButtonPushed());
            app.PlayPauseButton.Layout.Column = 1;

            app.StopButton = uibutton(transport, 'push', 'Text', 'Stop', 'Enable', 'off', ...
                'ButtonPushedFcn', @(~, ~) app.StopButtonPushed());
            app.StopButton.Layout.Column = 2;
        end

        function buildLatticeTab(app, parent)
            g = uigridlayout(parent, [1 2]);
            g.ColumnWidth = {'5x', '4x'};
            g.RowHeight = {'1x'};
            g.Padding = [8 8 8 8];
            g.ColumnSpacing = 8;

            app.PlotAx = uiaxes(g);
            app.PlotAx.Layout.Row = 1;
            app.PlotAx.Layout.Column = 1;
            title(app.PlotAx, 'Load a crystal or click "Build & Solve" to begin');
            app.PlotAx.DataAspectRatio = [1 1 1];
            grid(app.PlotAx, 'on');
            view(app.PlotAx, 3);

            app.SpecAx = uiaxes(g);
            app.SpecAx.Layout.Row = 1;
            app.SpecAx.Layout.Column = 2;
            title(app.SpecAx, 'Click "Build & Solve" to compute the spectrum');
            xlabel(app.SpecAx, 'Mode index');
            ylabel(app.SpecAx, '\omega (angular frequency)');
            grid(app.SpecAx, 'on');
        end

        function buildShearTab(app, parent)
            g = uigridlayout(parent, [2 3]);
            g.ColumnWidth = {340, '1x', '1x'};
            g.RowHeight = {'1x', '1x'};
            g.Padding = [8 8 8 8];
            g.RowSpacing = 8;
            g.ColumnSpacing = 8;

            panel = uipanel(g, 'Title', 'Shear Control Panel', 'BackgroundColor', 'w');
            panel.Layout.Row = [1 2];
            panel.Layout.Column = 1;
            app.buildShearControls(panel);

            app.UIAxes3D_NoForce = uiaxes(g);
            app.UIAxes3D_NoForce.Layout.Row = 1;
            app.UIAxes3D_NoForce.Layout.Column = 2;
            title(app.UIAxes3D_NoForce, 'Pristine Lattice (No Force)');

            app.UIAxes3D_WithForce = uiaxes(g);
            app.UIAxes3D_WithForce.Layout.Row = 1;
            app.UIAxes3D_WithForce.Layout.Column = 3;
            title(app.UIAxes3D_WithForce, 'Shear-Deformed Lattice (Force Distributed)');

            app.UIAxesMode_NoForce = uiaxes(g);
            app.UIAxesMode_NoForce.Layout.Row = 2;
            app.UIAxesMode_NoForce.Layout.Column = 2;
            title(app.UIAxesMode_NoForce, 'Normal Modes - Pristine');
            xlabel(app.UIAxesMode_NoForce, 'Mode Index');
            ylabel(app.UIAxesMode_NoForce, '\omega (angular frequency)');

            app.UIAxesMode_WithForce = uiaxes(g);
            app.UIAxesMode_WithForce.Layout.Row = 2;
            app.UIAxesMode_WithForce.Layout.Column = 3;
            title(app.UIAxesMode_WithForce, 'Normal Modes - Sheared');
            xlabel(app.UIAxesMode_WithForce, 'Mode Index');
            ylabel(app.UIAxesMode_WithForce, '\omega (angular frequency)');
        end

        function buildShearControls(app, parent)
            g = uigridlayout(parent, [12 1]);
            g.ColumnWidth = {'1x'};
            g.RowHeight = {22, 50, 34, 44, 26, 28, 34, 34, 34, 34, 96, '1x'};
            g.RowSpacing = 6;
            g.Padding = [10 10 10 10];
            g.Scrollable = 'on';

            lbl = uilabel(g, 'Text', 'Engineering Shear Strain (gamma):');
            lbl.Layout.Row = 1;

            app.ShearSlider = uislider(g, 'Limits', [0 0.5], 'Value', 0, ...
                'ValueChangedFcn', @(~, event) app.updateShearAnalysis(event.Value));
            app.ShearSlider.Layout.Row = 2;

            app.RelaxButton = uibutton(g, 'push', 'Text', 'Relax (Gauss-Jordan, Exp.5)', ...
                'FontWeight', 'bold', 'BackgroundColor', [0.20 0.45 0.80], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.RelaxButtonPushed());
            app.RelaxButton.Layout.Row = 3;

            app.RelaxStatusLabel = uilabel(g, ...
                'Text', 'Relax to solve for true equilibrium of interior atoms.', ...
                'FontColor', [0.4 0.4 0.4], 'WordWrap', 'on', 'VerticalAlignment', 'top');
            app.RelaxStatusLabel.Layout.Row = 4;

            hdr = uilabel(g, 'Text', 'Phonon softening & bond breaking (V7)', 'FontWeight', 'bold');
            hdr.Layout.Row = 5;
            hdr.Tooltip = 'Sweeps shear strain on the current lattice with three bond models (harmonic, +pre-stress, Morse with bond rupture).';

            gammaRow = uigridlayout(g, [1 2]);
            gammaRow.Layout.Row = 6;
            gammaRow.ColumnWidth = {'1x', 110};
            gammaRow.ColumnSpacing = 6;
            gammaRow.Padding = [0 0 0 0];
            lbl = uilabel(gammaRow, 'Text', 'Max shear strain gamma_max:');
            lbl.Layout.Column = 1;
            app.GammaMaxEdit = uieditfield(gammaRow, 'numeric', 'Limits', [0.05 1.5], 'Value', 0.5);
            app.GammaMaxEdit.Layout.Column = 2;

            app.PhononPlotsButton = uibutton(g, 'push', 'Text', 'Phonon Softening Plots (6 figures)', ...
                'FontWeight', 'bold', 'BackgroundColor', [0.47 0.67 0.19], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.runAnalysisPlots('phonon'));
            app.PhononPlotsButton.Layout.Row = 7;

            app.BondPlotsButton = uibutton(g, 'push', 'Text', 'Bond-Breaking Plots (2 figures)', ...
                'FontWeight', 'bold', 'BackgroundColor', [0.85 0.33 0.10], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.runAnalysisPlots('bond'));
            app.BondPlotsButton.Layout.Row = 8;

            app.DispersionButton = uibutton(g, 'push', 'Text', '3D Dispersion Surface (Brillouin zone)', ...
                'FontWeight', 'bold', 'BackgroundColor', [0.20 0.45 0.80], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.openDispersionSurface());
            app.DispersionButton.Layout.Row = 9;

            alphaRow = uigridlayout(g, [1 2]);
            alphaRow.Layout.Row = 10;
            alphaRow.ColumnWidth = {'1x', 80};
            alphaRow.ColumnSpacing = 6;
            alphaRow.Padding = [0 0 0 0];
            app.AlphaLabel = uilabel(alphaRow, 'Text', 'Morse alpha (critical stretch 17.3%):', 'WordWrap', 'on');
            app.AlphaLabel.Layout.Column = 1;
            app.AlphaEdit = uieditfield(alphaRow, 'numeric', 'Limits', [1 20], 'Value', 4, ...
                'Tooltip', 'ASSUMED default 4. Calibrate to the real material: critical stretch = ln2/alpha (ideal tensile strain).', ...
                'ValueChangedFcn', @(s, ~) set(app.AlphaLabel, 'Text', sprintf('Morse alpha (critical stretch %.1f%%):', 100 * log(2) / s.Value)));
            app.AlphaEdit.Layout.Column = 2;

            app.AnalysisStatusLabel = uilabel(g, ...
                'Text', ['Uses the lattice you built. The four axes on the right update live; the buttons above ' ...
                         'fill the Phonon Softening, Bond Breaking and Dispersion Surface tabs of this window.'], ...
                'WordWrap', 'on', 'VerticalAlignment', 'top', 'FontColor', [0.4 0.4 0.4]);
            app.AnalysisStatusLabel.Layout.Row = 11;
        end

        function resetAnalysisTabs(app)
            app.clearAnalysisTab(app.PhononTab, 'phonon');
            app.clearAnalysisTab(app.BondTab, 'bond');
            app.clearAnalysisTab(app.DispersionTab, 'dispersion');
        end

        function clearAnalysisTab(~, tab, kind)
            % Back to the "press the button" placeholder; also wipes a half-drawn tab after an error.
            switch kind
                case 'phonon'
                    msg = ['Press "Phonon Softening Plots" in the Shear Analysis tab.' newline newline ...
                           'The six plot groups P1-P6 (softening curve, DOS shift, spectral map, softest mode, ' ...
                           'stress vs phonon, transport) then appear here as sub-tabs.'];
                case 'bond'
                    msg = ['Press "Bond-Breaking Plots" in the Shear Analysis tab.' newline newline ...
                           'The two plot groups Q1 (where bonds break) and Q2 (when and why) then appear here as sub-tabs.'];
                otherwise
                    msg = ['Press "3D Dispersion Surface" in the Shear Analysis tab.' newline newline ...
                           'The Brillouin-zone surface and its own shear slider, bond-model menu and readout then appear here.'];
            end
            delete(allchild(tab));
            g = uigridlayout(tab, [1 1]);
            g.Padding = [30 30 30 30];
            lbl = uilabel(g, 'Text', msg, 'WordWrap', 'on', 'FontColor', [0.4 0.4 0.4], ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'center');
            lbl.Layout.Row = 1;
            lbl.Layout.Column = 1;
        end

        function setShearControlsEnabled(app, tf)
            if tf
                state = 'on';
            else
                state = 'off';
            end
            ctrls = {app.ShearSlider, app.RelaxButton, app.GammaMaxEdit, app.AlphaEdit, ...
                     app.PhononPlotsButton, app.BondPlotsButton, app.DispersionButton};
            for c = 1:numel(ctrls)
                if ~isempty(ctrls{c}) && isvalid(ctrls{c})
                    ctrls{c}.Enable = state;
                end
            end
        end

        function TabSelectionChanged(app, event)
            tab = event.NewValue;
            if tab == app.LatticeTab
                return;
            end
            if ~app.isBuilt
                app.TabGroup.SelectedTab = app.LatticeTab;
                app.showAlert('Please Build & Solve the lattice before running shear analysis.', 'Not Built');
                return;
            end
            if tab == app.ShearTab
                app.refreshShearViewIfNeeded();
            end
        end

        function refreshShearViewIfNeeded(app)
            % The four shear axes are only redrawn when the tab is actually looked at,
            % so Build & Solve stays as fast as it was when this lived in its own window.
            if ~app.shearViewDirty
                return;
            end
            app.shearViewDirty = false;
            app.updateShearAnalysis(app.ShearSlider.Value);
        end

        function openDispersionSurface(app)
            if ~app.isBuilt
                app.showAlert('Please Build & Solve the lattice first.', 'Not Built');
                return;
            end
            dlg = uiprogressdlg(app.UIFigure, 'Title', 'Dispersion surface', ...
                'Message', 'Sweeping the Brillouin zone under shear...', 'Indeterminate', 'on');
            cleanup = onCleanup(@() close(dlg)); %#ok<NASGU>
            try
                nb = size(app.AnalysisBasis, 1);
                brillouinShearSurface(app.AnalysisLV, app.AnalysisBasis, app.masses(1:nb), max(app.bonds(:, 3)), ...
                    app.BaseSpringConstant, app.GammaMaxEdit.Value, app.AlphaEdit.Value, app.DispersionTab);
                app.TabGroup.SelectedTab = app.DispersionTab;
                app.AnalysisStatusLabel.Text = 'Dispersion surface drawn in the "Dispersion Surface" tab; it has its own shear slider.';
            catch ME
                app.clearAnalysisTab(app.DispersionTab, 'dispersion');
                app.showAlert(ME.message, 'Dispersion surface error');
            end
        end

        function runAnalysisPlots(app, kind)
            if ~app.isBuilt
                app.showAlert('Please Build & Solve the lattice first.', 'Not Built');
                return;
            end
            dlg = uiprogressdlg(app.UIFigure, 'Title', 'Shear analysis', ...
                'Message', 'Sweeping shear strain...', 'Indeterminate', 'on');
            cleanup = onCleanup(@() close(dlg)); %#ok<NASGU>
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
                    plotShearSoftening(c, '', app.PhononTab);
                    app.TabGroup.SelectedTab = app.PhononTab;
                    msg = 'P1-P6 are in the "Phonon Softening" tab.';
                else
                    dlg.Message = 'Drawing bond-breaking plots...';
                    plotShearBonds(c, '', app.BondTab);
                    app.TabGroup.SelectedTab = app.BondTab;
                    msg = 'Q1-Q2 are in the "Bond Breaking" tab.';
                end
                R = c.cases{c.main};
                gb = min(R.gammaFirstBreak);
                if isnan(gb), gbs = 'none'; else, gbs = sprintf('%.3f', gb); end
                if isnan(R.gammaInstab), gi = 'none'; else, gi = sprintf('%.3f', R.gammaInstab); end
                app.AnalysisStatusLabel.Text = sprintf('%s First bond breaks at gamma = %s; lowest mode goes imaginary at gamma = %s.', msg, gbs, gi);
            catch ME
                if strcmp(kind, 'phonon')
                    app.clearAnalysisTab(app.PhononTab, 'phonon');
                else
                    app.clearAnalysisTab(app.BondTab, 'bond');
                end
                app.showAlert(ME.message, 'Analysis error');
            end
        end

        function updateShearAnalysis(app, shearVal)
            if isempty(app.pos) || isempty(app.BaseStiffnessMatrix)
                return;
            end

            [~, elementNames] = getElementDetails(app.masses);

            k = app.BaseSpringConstant;
            [~, ~, ~, freqs_base, ~] = applyShearForce(app.pos, app.bonds, app.masses, app.BaseStiffnessMatrix, k, 0.0);
            [coords_sh, force_dist, ~, freqs_sh, ~] = applyShearForce(app.pos, app.bonds, app.masses, app.BaseStiffnessMatrix, k, shearVal);

            cla(app.UIAxes3D_NoForce);
            title(app.UIAxes3D_NoForce, 'Pristine Lattice (No Force)');
            plotLattice(app.pos, app.bonds, app.UIAxes3D_NoForce);
            atomHover(app.UIAxes3D_NoForce, makeHoverLabels(elementNames, app.masses), @() app.pos);

            cla(app.UIAxes3D_WithForce);
            title(app.UIAxes3D_WithForce, 'Shear-Deformed Lattice (Force Distributed)');
            plotLatticeWithGradient(app.UIAxes3D_WithForce, coords_sh, elementNames, force_dist, app.bonds, 'Bond-force load (1 = load at gamma = 0.5)');
            atomHover(app.UIAxes3D_WithForce, makeHoverLabels(elementNames, app.masses), @() coords_sh);

            cla(app.UIAxesMode_NoForce);
            plot(app.UIAxesMode_NoForce, freqs_base, '-o', 'LineWidth', 1.5, 'Color', 'b');
            title(app.UIAxesMode_NoForce, 'Normal Modes - Pristine');
            xlabel(app.UIAxesMode_NoForce, 'Mode Index');
            ylabel(app.UIAxesMode_NoForce, '\omega (angular frequency)');

            yyaxis(app.UIAxesMode_WithForce, 'left');  cla(app.UIAxesMode_WithForce);
            yyaxis(app.UIAxesMode_WithForce, 'right');  cla(app.UIAxesMode_WithForce);
            yyaxis(app.UIAxesMode_WithForce, 'left');
            app.UIAxesMode_WithForce.YColor = [0 0 0];
            hold(app.UIAxesMode_WithForce, 'on');
            hPristine = plot(app.UIAxesMode_WithForce, freqs_base, '--', 'LineWidth', 1, 'Color', [0.6 0.6 0.6]);
            hSheared = plot(app.UIAxesMode_WithForce, freqs_sh, '-o', 'LineWidth', 1.5, 'Color', 'r');
            hold(app.UIAxesMode_WithForce, 'off');
            xlabel(app.UIAxesMode_WithForce, 'Mode Index');
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
            else
                title(app.UIAxesMode_WithForce, 'Normal Modes - Sheared');
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
            xlabel(app.UIAxesMode_WithForce, 'Mode Index');
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
                app.setStatus(sprintf('Loaded: %s (%d atoms/cell). Click "Build & Solve".', ...
                    file, size(app.basisFrac, 1)), [0.1 0.5 0.1]);
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

                app.setShearControlsEnabled(true);
                app.ShearSlider.Value = 0;
                app.RelaxStatusLabel.Text = 'Relax to solve for true equilibrium of interior atoms.';
                app.RelaxStatusLabel.FontColor = [0.4 0.4 0.4];
                app.shearViewDirty = true;
                app.resetAnalysisTabs();
                if app.TabGroup.SelectedTab == app.ShearTab
                    app.refreshShearViewIfNeeded();
                end

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
                app.PlayPauseButton.Enable = 'on';
                app.PlayPauseButton.Text = 'Play';
                app.PlayPauseButton.BackgroundColor = [0.30 0.70 0.35];
                app.StopButton.Enable = 'off';
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
            cla(app.SpecAx);
            hold(app.SpecAx, 'on');

            nModes = length(app.omega);
            idxAll = 1:nModes;
            plot(app.SpecAx, idxAll, app.omega, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor', [0.2 0.45 0.8]);
            hold(app.SpecAx, 'off');
            title(app.SpecAx, 'Frequency spectrum (red marker = selected mode)');
            xlabel(app.SpecAx, 'Mode index');
            ylabel(app.SpecAx, '\omega (angular frequency)');
            grid(app.SpecAx, 'on');
            app.SpecHighlightHandle = [];
        end

        function plotStaticLattice(app)
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

            ax.DataAspectRatio = [1 1 1];
            grid(ax, 'on');
            xlabel(ax, 'x'); ylabel(ax, 'y'); zlabel(ax, 'z');
            view(ax, 3);
            title(ax, sprintf('%d atoms, %d bonds (hover an atom for details)', size(app.pos, 1), nBonds));
            hold(ax, 'off');
        end

        function highlightModeOnSpectrum(app, modeIdx)
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

        function PlayPauseButtonPushed(app)
            if ~app.isBuilt
                return;
            end

            if ~isempty(app.animTimer) && isvalid(app.animTimer) && strcmp(app.animTimer.Running, 'on')
                stop(app.animTimer);
                app.PlayPauseButton.Text = 'Play';
                app.PlayPauseButton.BackgroundColor = [0.30 0.70 0.35];
                return;
            end

            if isempty(app.modeShape)
                app.ModeDropDownValueChanged();
            end

            app.TabGroup.SelectedTab = app.LatticeTab;

            if isempty(app.animTimer) || ~isvalid(app.animTimer)
                app.animTimer = timer( ...
                    'ExecutionMode', 'fixedRate', ...
                    'BusyMode', 'drop', ...
                    'Period', app.speedToPeriod(), ...
                    'TimerFcn', @(~, ~) app.animationTick());
            end
            start(app.animTimer);

            app.PlayPauseButton.Text = 'Pause';
            app.PlayPauseButton.BackgroundColor = [0.90 0.65 0.15];
            app.StopButton.Enable = 'on';
        end

        function StopButtonPushed(app)
            app.stopTimerIfRunning();
            app.animTime = 0;
            if app.isBuilt
                app.plotStaticLattice();
            end
            app.PlayPauseButton.Text = 'Play';
            app.PlayPauseButton.BackgroundColor = [0.30 0.70 0.35];
            app.PlayPauseButton.Enable = 'on';
            app.StopButton.Enable = 'off';
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
