classdef CrystalVibrationAppDesigner < matlab.apps.AppBase
    %CRYSTALVIBRATIONAPPDESIGNER App Designer version of CrystalVibrationApp.
    %   Same single-window GUI, written the way App Designer writes its code, so every
    %   callback below can be pasted straight into the matching callback of a .mlapp.
    %   See APPDESIGNER_NOTES.md next to this file for the build-it-in-Design-View recipe.
    %
    %   It also runs as-is from the command line:  CrystalVibrationAppDesigner
    %   It needs the other Version-8 .m files (autoRigidShells, buildDynamicalMatrix,
    %   jacobiEigenSolver, applyShearForce, relaxShearLattice, plotShearSoftening,
    %   plotShearBonds, brillouinShearSurface, shearPlotCanvas, ...) on the path.

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure              matlab.ui.Figure
        MainGrid              matlab.ui.container.GridLayout
        ControlPanel          matlab.ui.container.Panel
        ControlGrid           matlab.ui.container.GridLayout
        CrystalHeaderLabel    matlab.ui.control.Label
        LoadFileButton        matlab.ui.control.Button
        SupercellLabel        matlab.ui.control.Label
        NEditField            matlab.ui.control.NumericEditField
        kLabel                matlab.ui.control.Label
        kEditField            matlab.ui.control.NumericEditField
        BondNetworkLabel      matlab.ui.control.Label
        cutoffDropDown        matlab.ui.control.DropDown
        BuildButton           matlab.ui.control.Button
        AnimationHeaderLabel  matlab.ui.control.Label
        ModeLabel             matlab.ui.control.Label
        ModeDropDown          matlab.ui.control.DropDown
        AmplitudeLabel        matlab.ui.control.Label
        AmplitudeSlider       matlab.ui.control.Slider
        SpeedLabel            matlab.ui.control.Label
        SpeedSlider           matlab.ui.control.Slider
        TransportGrid         matlab.ui.container.GridLayout
        PlayPauseButton       matlab.ui.control.Button
        StopButton            matlab.ui.control.Button
        TabGroup              matlab.ui.container.TabGroup
        LatticeTab            matlab.ui.container.Tab
        LatticeGrid           matlab.ui.container.GridLayout
        PlotAx                matlab.ui.control.UIAxes
        SpecAx                matlab.ui.control.UIAxes
        ShearTab              matlab.ui.container.Tab
        ShearGrid             matlab.ui.container.GridLayout
        ShearControlPanel     matlab.ui.container.Panel
        ShearControlGrid      matlab.ui.container.GridLayout
        ShearStrainLabel      matlab.ui.control.Label
        ShearSlider           matlab.ui.control.Slider
        RelaxButton           matlab.ui.control.Button
        RelaxStatusLabel      matlab.ui.control.Label
        V7HeaderLabel         matlab.ui.control.Label
        GammaRowGrid          matlab.ui.container.GridLayout
        GammaMaxLabel         matlab.ui.control.Label
        GammaMaxEdit          matlab.ui.control.NumericEditField
        PhononPlotsButton     matlab.ui.control.Button
        BondPlotsButton       matlab.ui.control.Button
        DispersionButton      matlab.ui.control.Button
        AlphaRowGrid          matlab.ui.container.GridLayout
        AlphaLabel            matlab.ui.control.Label
        AlphaEdit             matlab.ui.control.NumericEditField
        AnalysisStatusLabel   matlab.ui.control.Label
        UIAxes3D_NoForce      matlab.ui.control.UIAxes
        UIAxes3D_WithForce    matlab.ui.control.UIAxes
        UIAxesMode_NoForce    matlab.ui.control.UIAxes
        UIAxesMode_WithForce  matlab.ui.control.UIAxes
        PhononTab             matlab.ui.container.Tab
        BondTab               matlab.ui.container.Tab
        DispersionTab         matlab.ui.container.Tab
        StatusLabel           matlab.ui.control.Label
    end

    % Model state. In App Designer: Code View -> "Property" -> Private property,
    % or just paste this whole block under the component properties block.
    properties (Access = private)
        latticeVectors              % lattice vectors of the loaded / built-in crystal
        basisFrac                   % fractional basis positions
        basisMasses                 % basis masses (amu)
        useCustomFile = false       % true once a CIF/txt file has been loaded

        pos                         % atom positions [Natoms x 3]
        bonds                       % bond list [i j d0 nx ny nz]
        masses                      % per-atom masses
        V                           % eigenvectors, columns sorted by omega
        omega                       % angular frequencies, ascending
        isBuilt = false             % true after a successful Build & Solve
        shearViewDirty = true       % shear axes need a redraw before being shown
        BaseStiffnessMatrix         % K of the unsheared lattice
        BaseSpringConstant          % k used for the current build

        modeShape                   % displacement pattern of the selected mode
        currentOmega = 0
        animTime = 0
        animTimer
        atomHandle                  % scatter3 of the animated atoms
        bondHandles                 % one line per bond
        SpecHighlightHandle         % red marker on the spectrum

        AnalysisLV                  % lattice vectors handed to the V7 sweeps
        AnalysisBasis               % basis handed to the V7 sweeps
        AnalysisCache               % last phononShearSuite result

        AnimFrameRateHz = 30        % animation redraw rate (Hz)
    end


    % Helper functions. In App Designer these go under Code View -> "Function" ->
    % Private function; the body of each is what you paste.
    methods (Access = private)

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

        function refreshShearViewIfNeeded(app)
            % The four shear axes are only redrawn when the tab is actually looked at,
            % so Build & Solve stays as fast as it was when this lived in its own window.
            if ~app.shearViewDirty
                return;
            end
            app.shearViewDirty = false;
            app.updateShearAnalysis(app.ShearSlider.Value);
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
            atomHover(app.UIAxes3D_NoForce, app.makeHoverLabels(elementNames, app.masses), @() app.pos);

            cla(app.UIAxes3D_WithForce);
            title(app.UIAxes3D_WithForce, 'Shear-Deformed Lattice (Force Distributed)');
            plotLatticeWithGradient(app.UIAxes3D_WithForce, coords_sh, elementNames, force_dist, app.bonds, 'Bond-force load (1 = load at gamma = 0.5)');
            atomHover(app.UIAxes3D_WithForce, app.makeHoverLabels(elementNames, app.masses), @() coords_sh);

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

            hoverTxt = app.makeHoverLabels(elementNames, app.masses);
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

        function labels = makeHoverLabels(~, elementNames, masses)
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

    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.resetAnalysisTabs();
            app.setShearControlsEnabled(false);
        end

        % Button pushed function: LoadFileButton
        function LoadFileButtonPushed(app, event)
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

        % Button pushed function: BuildButton
        function BuildButtonPushed(app, event)
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

        % Value changed function: ModeDropDown
        function ModeDropDownValueChanged(app, event)
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

        % Button pushed function: PlayPauseButton
        function PlayPauseButtonPushed(app, event)
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

        % Button pushed function: StopButton
        function StopButtonPushed(app, event)
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

        % Selection changed function: TabGroup
        function TabGroupSelectionChanged(app, event)
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

        % Value changed function: ShearSlider
        function ShearSliderValueChanged(app, event)
            app.updateShearAnalysis(app.ShearSlider.Value);
        end

        % Button pushed function: RelaxButton
        function RelaxButtonPushed(app, event)
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
            atomHover(app.UIAxes3D_WithForce, app.makeHoverLabels(elementNames, app.masses), @() coords_relaxed);
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

        % Button pushed function: PhononPlotsButton
        function PhononPlotsButtonPushed(app, event)
            app.runAnalysisPlots('phonon');
        end

        % Button pushed function: BondPlotsButton
        function BondPlotsButtonPushed(app, event)
            app.runAnalysisPlots('bond');
        end

        % Button pushed function: DispersionButton
        function DispersionButtonPushed(app, event)
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

        % Value changed function: AlphaEdit
        function AlphaEditValueChanged(app, event)
            app.AlphaLabel.Text = sprintf('Morse alpha (critical stretch %.1f%%):', 100 * log(2) / app.AlphaEdit.Value);
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            app.stopTimerIfRunning();
            delete(app)
        end

    end

    % Component initialization.
    % App Designer GENERATES this section and locks it in Code View - you cannot paste it there.
    % Build the same tree in Design View instead (see APPDESIGNER_NOTES.md); this is the spec
    % of exactly which component goes where, with which properties.
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [60 60 1500 880];
            app.UIFigure.Name = 'Crystal Vibration Explorer';
            app.UIFigure.Color = [1 1 1];
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create MainGrid
            app.MainGrid = uigridlayout(app.UIFigure);
            app.MainGrid.ColumnWidth = {340, '1x'};
            app.MainGrid.RowHeight = {'1x', 24};
            app.MainGrid.RowSpacing = 6;
            app.MainGrid.ColumnSpacing = 8;

            % Create ControlPanel
            app.ControlPanel = uipanel(app.MainGrid);
            app.ControlPanel.Title = 'Model & Animation';
            app.ControlPanel.BackgroundColor = [1 1 1];
            app.ControlPanel.Layout.Row = 1;
            app.ControlPanel.Layout.Column = 1;

            % Create ControlGrid
            app.ControlGrid = uigridlayout(app.ControlPanel);
            app.ControlGrid.ColumnWidth = {'1x'};
            app.ControlGrid.RowHeight = {24, 32, 22, 24, 22, 24, 22, 24, 38, 26, 22, 24, 22, 50, 22, 50, 34, '1x'};
            app.ControlGrid.RowSpacing = 6;
            app.ControlGrid.Padding = [10 10 10 10];
            app.ControlGrid.Scrollable = 'on';

            % Create CrystalHeaderLabel
            app.CrystalHeaderLabel = uilabel(app.ControlGrid);
            app.CrystalHeaderLabel.FontWeight = 'bold';
            app.CrystalHeaderLabel.Layout.Row = 1;
            app.CrystalHeaderLabel.Text = 'Crystal & bond network';

            % Create LoadFileButton
            app.LoadFileButton = uibutton(app.ControlGrid, 'push');
            app.LoadFileButton.ButtonPushedFcn = createCallbackFcn(app, @LoadFileButtonPushed, true);
            app.LoadFileButton.Layout.Row = 2;
            app.LoadFileButton.Text = 'Load Crystal File...';

            % Create SupercellLabel
            app.SupercellLabel = uilabel(app.ControlGrid);
            app.SupercellLabel.Layout.Row = 3;
            app.SupercellLabel.Text = 'Supercell size (fixed): 1 x 1 x 1 cell';

            % Create NEditField
            app.NEditField = uieditfield(app.ControlGrid, 'numeric');
            app.NEditField.Limits = [1 1];
            app.NEditField.Editable = 'off';
            app.NEditField.Layout.Row = 4;
            app.NEditField.Value = 1;

            % Create kLabel
            app.kLabel = uilabel(app.ControlGrid);
            app.kLabel.Layout.Row = 5;
            app.kLabel.Text = 'Spring constant k:';

            % Create kEditField
            % k = 1.0 is a reduced-unit default: use a force constant fitted to the real material.
            app.kEditField = uieditfield(app.ControlGrid, 'numeric');
            app.kEditField.Limits = [0.001 Inf];
            app.kEditField.Layout.Row = 6;
            app.kEditField.Value = 1;

            % Create BondNetworkLabel
            app.BondNetworkLabel = uilabel(app.ControlGrid);
            app.BondNetworkLabel.Layout.Row = 7;
            app.BondNetworkLabel.Text = 'Bond network:';

            % Create cutoffDropDown
            app.cutoffDropDown = uidropdown(app.ControlGrid);
            app.cutoffDropDown.Items = {'Nearest neighbors only', '1st + 2nd neighbors (fixed)', 'Auto (grows shells until rigid)'};
            app.cutoffDropDown.Layout.Row = 8;
            app.cutoffDropDown.Value = 'Auto (grows shells until rigid)';

            % Create BuildButton
            app.BuildButton = uibutton(app.ControlGrid, 'push');
            app.BuildButton.ButtonPushedFcn = createCallbackFcn(app, @BuildButtonPushed, true);
            app.BuildButton.BackgroundColor = [0.2 0.45 0.8];
            app.BuildButton.FontWeight = 'bold';
            app.BuildButton.FontColor = [1 1 1];
            app.BuildButton.Layout.Row = 9;
            app.BuildButton.Text = 'Build && Solve';

            % Create AnimationHeaderLabel
            app.AnimationHeaderLabel = uilabel(app.ControlGrid);
            app.AnimationHeaderLabel.FontWeight = 'bold';
            app.AnimationHeaderLabel.Layout.Row = 10;
            app.AnimationHeaderLabel.Text = 'Mode animation';

            % Create ModeLabel
            app.ModeLabel = uilabel(app.ControlGrid);
            app.ModeLabel.Layout.Row = 11;
            app.ModeLabel.Text = 'Vibrational mode:';

            % Create ModeDropDown
            app.ModeDropDown = uidropdown(app.ControlGrid);
            app.ModeDropDown.Items = {'(build the model first)'};
            app.ModeDropDown.ValueChangedFcn = createCallbackFcn(app, @ModeDropDownValueChanged, true);
            app.ModeDropDown.Enable = 'off';
            app.ModeDropDown.Layout.Row = 12;
            app.ModeDropDown.Value = '(build the model first)';

            % Create AmplitudeLabel
            app.AmplitudeLabel = uilabel(app.ControlGrid);
            app.AmplitudeLabel.Layout.Row = 13;
            app.AmplitudeLabel.Text = 'Amplitude:';

            % Create AmplitudeSlider
            app.AmplitudeSlider = uislider(app.ControlGrid);
            app.AmplitudeSlider.Limits = [0 1];
            app.AmplitudeSlider.Layout.Row = 14;
            app.AmplitudeSlider.Value = 0.25;

            % Create SpeedLabel
            app.SpeedLabel = uilabel(app.ControlGrid);
            app.SpeedLabel.Layout.Row = 15;
            app.SpeedLabel.Text = 'Playback speed (cycles/sec):';

            % Create SpeedSlider
            app.SpeedSlider = uislider(app.ControlGrid);
            app.SpeedSlider.Limits = [0.1 10];
            app.SpeedSlider.Layout.Row = 16;
            app.SpeedSlider.Value = 2;

            % Create TransportGrid
            app.TransportGrid = uigridlayout(app.ControlGrid);
            app.TransportGrid.ColumnWidth = {'1x', '1x'};
            app.TransportGrid.RowHeight = {'1x'};
            app.TransportGrid.ColumnSpacing = 8;
            app.TransportGrid.Padding = [0 0 0 0];
            app.TransportGrid.Layout.Row = 17;

            % Create PlayPauseButton
            app.PlayPauseButton = uibutton(app.TransportGrid, 'push');
            app.PlayPauseButton.ButtonPushedFcn = createCallbackFcn(app, @PlayPauseButtonPushed, true);
            app.PlayPauseButton.BackgroundColor = [0.3 0.7 0.35];
            app.PlayPauseButton.FontColor = [1 1 1];
            app.PlayPauseButton.Enable = 'off';
            app.PlayPauseButton.Layout.Column = 1;
            app.PlayPauseButton.Text = 'Play';

            % Create StopButton
            app.StopButton = uibutton(app.TransportGrid, 'push');
            app.StopButton.ButtonPushedFcn = createCallbackFcn(app, @StopButtonPushed, true);
            app.StopButton.Enable = 'off';
            app.StopButton.Layout.Column = 2;
            app.StopButton.Text = 'Stop';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.MainGrid);
            app.TabGroup.SelectionChangedFcn = createCallbackFcn(app, @TabGroupSelectionChanged, true);
            app.TabGroup.Layout.Row = 1;
            app.TabGroup.Layout.Column = 2;

            % Create LatticeTab
            app.LatticeTab = uitab(app.TabGroup);
            app.LatticeTab.Title = 'Lattice & Modes';

            % Create LatticeGrid
            app.LatticeGrid = uigridlayout(app.LatticeTab);
            app.LatticeGrid.ColumnWidth = {'5x', '4x'};
            app.LatticeGrid.RowHeight = {'1x'};
            app.LatticeGrid.ColumnSpacing = 8;
            app.LatticeGrid.Padding = [8 8 8 8];

            % Create PlotAx
            app.PlotAx = uiaxes(app.LatticeGrid);
            title(app.PlotAx, 'Load a crystal or click "Build & Solve" to begin')
            xlabel(app.PlotAx, 'x')
            ylabel(app.PlotAx, 'y')
            zlabel(app.PlotAx, 'z')
            app.PlotAx.DataAspectRatio = [1 1 1];
            app.PlotAx.View = [-37.5 30];
            app.PlotAx.XGrid = 'on';
            app.PlotAx.YGrid = 'on';
            app.PlotAx.ZGrid = 'on';
            app.PlotAx.Layout.Row = 1;
            app.PlotAx.Layout.Column = 1;

            % Create SpecAx
            app.SpecAx = uiaxes(app.LatticeGrid);
            title(app.SpecAx, 'Click "Build & Solve" to compute the spectrum')
            xlabel(app.SpecAx, 'Mode index')
            ylabel(app.SpecAx, '\omega (angular frequency)')
            app.SpecAx.XGrid = 'on';
            app.SpecAx.YGrid = 'on';
            app.SpecAx.Layout.Row = 1;
            app.SpecAx.Layout.Column = 2;

            % Create ShearTab
            app.ShearTab = uitab(app.TabGroup);
            app.ShearTab.Title = 'Shear Analysis';

            % Create ShearGrid
            app.ShearGrid = uigridlayout(app.ShearTab);
            app.ShearGrid.ColumnWidth = {340, '1x', '1x'};
            app.ShearGrid.RowHeight = {'1x', '1x'};
            app.ShearGrid.RowSpacing = 8;
            app.ShearGrid.ColumnSpacing = 8;
            app.ShearGrid.Padding = [8 8 8 8];

            % Create ShearControlPanel
            app.ShearControlPanel = uipanel(app.ShearGrid);
            app.ShearControlPanel.Title = 'Shear Control Panel';
            app.ShearControlPanel.BackgroundColor = [1 1 1];
            app.ShearControlPanel.Layout.Row = [1 2];
            app.ShearControlPanel.Layout.Column = 1;

            % Create ShearControlGrid
            app.ShearControlGrid = uigridlayout(app.ShearControlPanel);
            app.ShearControlGrid.ColumnWidth = {'1x'};
            app.ShearControlGrid.RowHeight = {22, 50, 34, 44, 26, 28, 34, 34, 34, 34, 96, '1x'};
            app.ShearControlGrid.RowSpacing = 6;
            app.ShearControlGrid.Padding = [10 10 10 10];
            app.ShearControlGrid.Scrollable = 'on';

            % Create ShearStrainLabel
            app.ShearStrainLabel = uilabel(app.ShearControlGrid);
            app.ShearStrainLabel.Layout.Row = 1;
            app.ShearStrainLabel.Text = 'Engineering Shear Strain (gamma):';

            % Create ShearSlider
            app.ShearSlider = uislider(app.ShearControlGrid);
            app.ShearSlider.Limits = [0 0.5];
            app.ShearSlider.ValueChangedFcn = createCallbackFcn(app, @ShearSliderValueChanged, true);
            app.ShearSlider.Layout.Row = 2;
            app.ShearSlider.Value = 0;

            % Create RelaxButton
            app.RelaxButton = uibutton(app.ShearControlGrid, 'push');
            app.RelaxButton.ButtonPushedFcn = createCallbackFcn(app, @RelaxButtonPushed, true);
            app.RelaxButton.BackgroundColor = [0.2 0.45 0.8];
            app.RelaxButton.FontWeight = 'bold';
            app.RelaxButton.FontColor = [1 1 1];
            app.RelaxButton.Layout.Row = 3;
            app.RelaxButton.Text = 'Relax (Gauss-Jordan, Exp.5)';

            % Create RelaxStatusLabel
            app.RelaxStatusLabel = uilabel(app.ShearControlGrid);
            app.RelaxStatusLabel.VerticalAlignment = 'top';
            app.RelaxStatusLabel.WordWrap = 'on';
            app.RelaxStatusLabel.FontColor = [0.4 0.4 0.4];
            app.RelaxStatusLabel.Layout.Row = 4;
            app.RelaxStatusLabel.Text = 'Relax to solve for true equilibrium of interior atoms.';

            % Create V7HeaderLabel
            app.V7HeaderLabel = uilabel(app.ShearControlGrid);
            app.V7HeaderLabel.FontWeight = 'bold';
            app.V7HeaderLabel.Tooltip = {'Sweeps shear strain on the current lattice with three bond models (harmonic, +pre-stress, Morse with bond rupture).'};
            app.V7HeaderLabel.Layout.Row = 5;
            app.V7HeaderLabel.Text = 'Phonon softening & bond breaking (V7)';

            % Create GammaRowGrid
            app.GammaRowGrid = uigridlayout(app.ShearControlGrid);
            app.GammaRowGrid.ColumnWidth = {'1x', 110};
            app.GammaRowGrid.RowHeight = {'1x'};
            app.GammaRowGrid.ColumnSpacing = 6;
            app.GammaRowGrid.Padding = [0 0 0 0];
            app.GammaRowGrid.Layout.Row = 6;

            % Create GammaMaxLabel
            app.GammaMaxLabel = uilabel(app.GammaRowGrid);
            app.GammaMaxLabel.Layout.Column = 1;
            app.GammaMaxLabel.Text = 'Max shear strain gamma_max:';

            % Create GammaMaxEdit
            app.GammaMaxEdit = uieditfield(app.GammaRowGrid, 'numeric');
            app.GammaMaxEdit.Limits = [0.05 1.5];
            app.GammaMaxEdit.Layout.Column = 2;
            app.GammaMaxEdit.Value = 0.5;

            % Create PhononPlotsButton
            app.PhononPlotsButton = uibutton(app.ShearControlGrid, 'push');
            app.PhononPlotsButton.ButtonPushedFcn = createCallbackFcn(app, @PhononPlotsButtonPushed, true);
            app.PhononPlotsButton.BackgroundColor = [0.47 0.67 0.19];
            app.PhononPlotsButton.FontWeight = 'bold';
            app.PhononPlotsButton.FontColor = [1 1 1];
            app.PhononPlotsButton.Layout.Row = 7;
            app.PhononPlotsButton.Text = 'Phonon Softening Plots (6 figures)';

            % Create BondPlotsButton
            app.BondPlotsButton = uibutton(app.ShearControlGrid, 'push');
            app.BondPlotsButton.ButtonPushedFcn = createCallbackFcn(app, @BondPlotsButtonPushed, true);
            app.BondPlotsButton.BackgroundColor = [0.85 0.33 0.1];
            app.BondPlotsButton.FontWeight = 'bold';
            app.BondPlotsButton.FontColor = [1 1 1];
            app.BondPlotsButton.Layout.Row = 8;
            app.BondPlotsButton.Text = 'Bond-Breaking Plots (2 figures)';

            % Create DispersionButton
            app.DispersionButton = uibutton(app.ShearControlGrid, 'push');
            app.DispersionButton.ButtonPushedFcn = createCallbackFcn(app, @DispersionButtonPushed, true);
            app.DispersionButton.BackgroundColor = [0.2 0.45 0.8];
            app.DispersionButton.FontWeight = 'bold';
            app.DispersionButton.FontColor = [1 1 1];
            app.DispersionButton.Layout.Row = 9;
            app.DispersionButton.Text = '3D Dispersion Surface (Brillouin zone)';

            % Create AlphaRowGrid
            app.AlphaRowGrid = uigridlayout(app.ShearControlGrid);
            app.AlphaRowGrid.ColumnWidth = {'1x', 80};
            app.AlphaRowGrid.RowHeight = {'1x'};
            app.AlphaRowGrid.ColumnSpacing = 6;
            app.AlphaRowGrid.Padding = [0 0 0 0];
            app.AlphaRowGrid.Layout.Row = 10;

            % Create AlphaLabel
            app.AlphaLabel = uilabel(app.AlphaRowGrid);
            app.AlphaLabel.WordWrap = 'on';
            app.AlphaLabel.Layout.Column = 1;
            app.AlphaLabel.Text = 'Morse alpha (critical stretch 17.3%):';

            % Create AlphaEdit
            app.AlphaEdit = uieditfield(app.AlphaRowGrid, 'numeric');
            app.AlphaEdit.Limits = [1 20];
            app.AlphaEdit.ValueChangedFcn = createCallbackFcn(app, @AlphaEditValueChanged, true);
            app.AlphaEdit.Tooltip = {'ASSUMED default 4. Calibrate to the real material: critical stretch = ln2/alpha (ideal tensile strain).'};
            app.AlphaEdit.Layout.Column = 2;
            app.AlphaEdit.Value = 4;

            % Create AnalysisStatusLabel
            app.AnalysisStatusLabel = uilabel(app.ShearControlGrid);
            app.AnalysisStatusLabel.VerticalAlignment = 'top';
            app.AnalysisStatusLabel.WordWrap = 'on';
            app.AnalysisStatusLabel.FontColor = [0.4 0.4 0.4];
            app.AnalysisStatusLabel.Layout.Row = 11;
            app.AnalysisStatusLabel.Text = 'Uses the lattice you built. The four axes on the right update live; the buttons above fill the Phonon Softening, Bond Breaking and Dispersion Surface tabs of this window.';

            % Create UIAxes3D_NoForce
            app.UIAxes3D_NoForce = uiaxes(app.ShearGrid);
            title(app.UIAxes3D_NoForce, 'Pristine Lattice (No Force)')
            app.UIAxes3D_NoForce.Layout.Row = 1;
            app.UIAxes3D_NoForce.Layout.Column = 2;

            % Create UIAxes3D_WithForce
            app.UIAxes3D_WithForce = uiaxes(app.ShearGrid);
            title(app.UIAxes3D_WithForce, 'Shear-Deformed Lattice (Force Distributed)')
            app.UIAxes3D_WithForce.Layout.Row = 1;
            app.UIAxes3D_WithForce.Layout.Column = 3;

            % Create UIAxesMode_NoForce
            app.UIAxesMode_NoForce = uiaxes(app.ShearGrid);
            title(app.UIAxesMode_NoForce, 'Normal Modes - Pristine')
            xlabel(app.UIAxesMode_NoForce, 'Mode Index')
            ylabel(app.UIAxesMode_NoForce, '\omega (angular frequency)')
            app.UIAxesMode_NoForce.Layout.Row = 2;
            app.UIAxesMode_NoForce.Layout.Column = 2;

            % Create UIAxesMode_WithForce
            app.UIAxesMode_WithForce = uiaxes(app.ShearGrid);
            title(app.UIAxesMode_WithForce, 'Normal Modes - Sheared')
            xlabel(app.UIAxesMode_WithForce, 'Mode Index')
            ylabel(app.UIAxesMode_WithForce, '\omega (angular frequency)')
            app.UIAxesMode_WithForce.Layout.Row = 2;
            app.UIAxesMode_WithForce.Layout.Column = 3;

            % Create PhononTab
            app.PhononTab = uitab(app.TabGroup);
            app.PhononTab.Title = 'Phonon Softening';
            app.PhononTab.BackgroundColor = [1 1 1];

            % Create BondTab
            app.BondTab = uitab(app.TabGroup);
            app.BondTab.Title = 'Bond Breaking';
            app.BondTab.BackgroundColor = [1 1 1];

            % Create DispersionTab
            app.DispersionTab = uitab(app.TabGroup);
            app.DispersionTab.Title = 'Dispersion Surface';
            app.DispersionTab.BackgroundColor = [1 1 1];

            % Create StatusLabel
            app.StatusLabel = uilabel(app.MainGrid);
            app.StatusLabel.FontColor = [0.4 0.4 0.4];
            app.StatusLabel.Layout.Row = 2;
            app.StatusLabel.Layout.Column = [1 2];
            app.StatusLabel.Text = 'No file loaded (using default bcc structure).';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = CrystalVibrationAppDesigner

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
