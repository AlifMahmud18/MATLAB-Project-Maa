classdef CrystalVibrationApp < handle
    %CRYSTALVIBRATIONAPP Interactive crystal vibration viewer with element colors and hover tooltips.

    properties (Access = private)
        % --- UI components ---
        UIFigure
        LoadFileButton
        StatusLabel
        StructureDropDown
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

        % --- Combined plot window (3D view + frequency spectrum) ---
        PlotFig
        Panel3D
        PanelSpec
        PlotAx
        SpecAx
        View3DButton
        ViewSpecButton
        SpecHighlightHandle

        % --- Crystal definition currently loaded ---
        latticeVectors
        basisFrac
        basisMasses
        useCustomFile = false

        % --- Solved model ---
        pos
        bonds
        masses         % array of masses for all supercell atoms
        V              % eigenvector matrix
        omega          % vector of angular frequencies
        isBuilt = false

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
                'FontSize', 16, 'FontWeight', 'bold', 'Position', [20 615 280 25]);

            % --- Load / structure selection ---
            app.LoadFileButton = uibutton(app.UIFigure, 'push', ...
                'Text', 'Load Crystal File...', 'Position', [20 580 260 30], ...
                'ButtonPushedFcn', @(~, ~) app.LoadFileButtonPushed());

            app.StatusLabel = uilabel(app.UIFigure, 'Text', 'No file loaded (using built-in structure)', ...
                'Position', [20 555 280 20], 'FontColor', [0.4 0.4 0.4]);

            uilabel(app.UIFigure, 'Text', 'Built-in structure:', 'Position', [20 525 260 20]);
            app.StructureDropDown = uidropdown(app.UIFigure, ...
                'Items', {'sc (simple cubic)', 'bcc', 'fcc'}, 'Value', 'bcc', ...
                'Position', [20 500 260 25], ...
                'ValueChangedFcn', @(~, ~) app.StructureDropDownValueChanged());

            % --- Physics parameters ---
            uilabel(app.UIFigure, 'Text', 'Supercell size N (NxNxN cells):', 'Position', [20 465 260 20]);
            app.NEditField = uieditfield(app.UIFigure, 'numeric', ...
                'Value', 2, 'Limits', [1 4], ...
                'Position', [20 440 260 25]);

            uilabel(app.UIFigure, 'Text', 'Spring constant k:', 'Position', [20 410 260 20]);
            app.kEditField = uieditfield(app.UIFigure, 'numeric', ...
                'Value', 1.0, 'Limits', [0.001 Inf], 'Position', [20 385 260 25]);

            uilabel(app.UIFigure, 'Text', 'Bond network:', 'Position', [20 350 260 20]);
            app.cutoffDropDown = uidropdown(app.UIFigure, ...
                'Items', {'Nearest neighbors only', '1st + 2nd neighbors (recommended, rigid)'}, ...
                'Value', '1st + 2nd neighbors (recommended, rigid)', ...
                'Position', [20 325 260 25]);

            app.BuildButton = uibutton(app.UIFigure, 'push', ...
                'Text', 'Build && Solve', 'Position', [20 280 260 35], ...
                'FontWeight', 'bold', 'BackgroundColor', [0.20 0.45 0.80], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.BuildButtonPushed());

            % --- Mode selection ---
            uilabel(app.UIFigure, 'Text', 'Vibrational mode:', 'Position', [20 245 260 20]);
            app.ModeDropDown = uidropdown(app.UIFigure, ...
                'Items', {'(build the model first)'}, 'Enable', 'off', ...
                'Position', [20 220 260 25], ...
                'ValueChangedFcn', @(~, ~) app.ModeDropDownValueChanged());

            % --- Animation controls ---
            uilabel(app.UIFigure, 'Text', 'Amplitude:', 'Position', [20 185 260 20]);
            app.AmplitudeSlider = uislider(app.UIFigure, ...
                'Limits', [0 1], 'Value', 0.25, 'Position', [25 170 250 3]);

            uilabel(app.UIFigure, 'Text', 'Speed (fps):', 'Position', [20 125 260 20]);
            app.SpeedSlider = uislider(app.UIFigure, ...
                'Limits', [5 60], 'Value', 30, 'Position', [25 110 250 3], ...
                'ValueChangedFcn', @(~, ~) app.SpeedSliderValueChanged());

            app.PlayButton = uibutton(app.UIFigure, 'push', 'Text', 'Play', ...
                'Position', [20 55 80 32], 'Enable', 'off', ...
                'BackgroundColor', [0.30 0.70 0.35], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.PlayButtonPushed());
            app.PauseButton = uibutton(app.UIFigure, 'push', 'Text', 'Pause', ...
                'Position', [110 55 80 32], 'Enable', 'off', ...
                'ButtonPushedFcn', @(~, ~) app.PauseButtonPushed());
            app.StopButton = uibutton(app.UIFigure, 'push', 'Text', 'Stop', ...
                'Position', [200 55 80 32], 'Enable', 'off', ...
                'ButtonPushedFcn', @(~, ~) app.StopButtonPushed());

            uilabel(app.UIFigure, 'Text', 'The 3D view and frequency spectrum open together in a separate window.', ...
                'Position', [20 15 280 30], 'FontColor', [0.4 0.4 0.4]);

            app.ensurePlotFigure();
        end

        function ensurePlotFigure(app)
            if ~isempty(app.PlotFig) && isvalid(app.PlotFig)
                return;
            end

            app.PlotFig = figure('Name', 'Crystal Vibration Explorer - Views', ...
                'NumberTitle', 'off', 'Position', [440 80 680 700], 'Color', 'w');

            % --- Toggle buttons to switch between the two views ---
            app.View3DButton = uicontrol(app.PlotFig, 'Style', 'pushbutton', ...
                'String', '3D View', 'Position', [10 654 200 34], ...
                'FontWeight', 'bold', 'Callback', @(~, ~) app.switchToView('3D'));
            app.ViewSpecButton = uicontrol(app.PlotFig, 'Style', 'pushbutton', ...
                'String', 'Frequency Spectrum', 'Position', [220 654 200 34], ...
                'FontWeight', 'bold', 'Callback', @(~, ~) app.switchToView('Spec'));

            % --- Panel for the 3D view ---
            app.Panel3D = uipanel(app.PlotFig, 'Units', 'pixels', ...
                'Position', [10 10 660 630], 'BorderType', 'none', 'BackgroundColor', 'w');
            app.PlotAx = axes('Parent', app.Panel3D);
            title(app.PlotAx, 'Load a crystal and click "Build & Solve" to begin');
            axis(app.PlotAx, 'equal');
            grid(app.PlotAx, 'on');
            view(app.PlotAx, 3);

            % --- Panel for the frequency spectrum ---
            app.PanelSpec = uipanel(app.PlotFig, 'Units', 'pixels', ...
                'Position', [10 10 660 630], 'BorderType', 'none', 'BackgroundColor', 'w');
            app.SpecAx = axes('Parent', app.PanelSpec);
            title(app.SpecAx, 'Click "Build & Solve" to compute the spectrum');
            xlabel(app.SpecAx, 'Mode index');
            ylabel(app.SpecAx, '\omega (angular frequency)');
            grid(app.SpecAx, 'on');

            app.switchToView('3D');
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

        function StructureDropDownValueChanged(app)
            app.useCustomFile = false;
            app.StatusLabel.Text = 'Using built-in structure (no file loaded)';
            app.StatusLabel.FontColor = [0.4 0.4 0.4];
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
                    switch app.StructureDropDown.Value
                        case 'sc (simple cubic)'
                            [latVec, basis, bMasses] = getBuiltinCrystalDef('sc');
                        case 'bcc'
                            [latVec, basis, bMasses] = getBuiltinCrystalDef('bcc');
                        case 'fcc'
                            [latVec, basis, bMasses] = getBuiltinCrystalDef('fcc');
                    end
                end

                app.setStatus('Building lattice...', [0.2 0.2 0.6]);
                [app.pos, app.bonds, app.masses] = generateLatticeGeneral(latVec, basis, bMasses, N, numShells);

                app.setStatus('Assembling dynamical matrix...', [0.2 0.2 0.6]);
                [~, ~, Dmat] = buildDynamicalMatrix(app.pos, app.bonds, k, app.masses);

                app.setStatus('Solving eigenproblem...', [0.2 0.2 0.6]);
                [V, eigVals] = jacobiEigenSolver(Dmat, 1e-12, 200);
                [eigVals, order] = sort(eigVals);
                app.V = V(:, order);
                app.omega = sqrt(max(eigVals, 0));

                app.populateModeDropdown();
                app.plotStaticLattice();
                app.plotFrequencySpectrum();

                Natoms = size(app.pos, 1);
                numZero = sum(eigVals < 1e-6);
                app.setStatus(sprintf('Done: %d atoms, %d bonds, %d modes (%d rigid-body zero modes).', ...
                    Natoms, size(app.bonds, 1), length(app.omega), numZero), [0.1 0.5 0.1]);

                app.isBuilt = true;
                app.ModeDropDown.Enable = 'on';
                app.PlayButton.Enable = 'on';
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
            
            % Set items and string value directly without using ItemsData
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
            isZero = app.omega < 1e-6;

            % Rigid-body / zero modes in gray, vibrating modes in blue
            stem(app.SpecAx, idxAll(isZero), app.omega(isZero), ...
                'Color', [0.6 0.6 0.6], 'Marker', 'o', 'MarkerFaceColor', [0.6 0.6 0.6], ...
                'LineWidth', 1.0, 'DisplayName', 'Rigid-body (\omega \approx 0)');
            stem(app.SpecAx, idxAll(~isZero), app.omega(~isZero), ...
                'Color', [0.20 0.45 0.80], 'Marker', 'o', 'MarkerFaceColor', [0.20 0.45 0.80], ...
                'LineWidth', 1.0, 'DisplayName', 'Vibrational modes');

            % Highlight the currently selected mode, if any
            app.SpecHighlightHandle = plot(app.SpecAx, nan, nan, 'o', ...
                'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'LineWidth', 2, ...
                'DisplayName', 'Selected mode');

            xlabel(app.SpecAx, 'Mode index');
            ylabel(app.SpecAx, '\omega (angular frequency)');
            title(app.SpecAx, sprintf('Frequency spectrum (%d modes, %d rigid-body)', ...
                nModes, sum(isZero)));
            legend(app.SpecAx, 'Location', 'northwest');
            grid(app.SpecAx, 'on');
            xlim(app.SpecAx, [0.5, nModes + 0.5]);
            hold(app.SpecAx, 'off');

            app.updateSpectrumHighlight();
        end

        function updateSpectrumHighlight(app)
            if isempty(app.SpecHighlightHandle) || ~isvalid(app.SpecHighlightHandle)
                return;
            end
            if isempty(app.ModeDropDown.Value) || ~app.isBuilt
                set(app.SpecHighlightHandle, 'XData', nan, 'YData', nan);
                return;
            end
            modeIdx = app.getSelectedModeIndex();
            set(app.SpecHighlightHandle, 'XData', modeIdx, 'YData', app.omega(modeIdx));
        end

        function ModeDropDownValueChanged(app)
            app.updateSpectrumHighlight();
        end

        function plotStaticLattice(app)
            app.stopTimerIfRunning();
            app.ensurePlotFigure();
            cla(app.PlotAx);
            hold(app.PlotAx, 'on');

            % 1. Draw bonds
            Nbonds = size(app.bonds, 1);
            for kk = 1:Nbonds
                i = app.bonds(kk, 1);
                j = app.bonds(kk, 2);
                plot3(app.PlotAx, ...
                    [app.pos(i,1) app.pos(j,1)], [app.pos(i,2) app.pos(j,2)], [app.pos(i,3) app.pos(j,3)], ...
                    'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
            end

            % 2. Map elements to colors and names
            [atomColors, elementNames] = getElementDetails(app.masses);

            % 3. Plot colored atoms
            app.atomHandle = scatter3(app.PlotAx, app.pos(:,1), app.pos(:,2), app.pos(:,3), ...
                160, atomColors, 'filled', 'MarkerEdgeColor', 'k');

            % 4. Add hover / click data tips safely (R2019b+)
            if isprop(app.atomHandle, 'DataTipTemplate') && exist('dataTipTextRow', 'file')
                dt = app.atomHandle.DataTipTemplate;
                dt.DataTipRows = [ ...
                    dataTipTextRow('Element', elementNames); ...
                    dataTipTextRow('Mass', app.masses, '%.3f amu'); ...
                    dataTipTextRow('X', app.pos(:,1), '%.3f'); ...
                    dataTipTextRow('Y', app.pos(:,2), '%.3f'); ...
                    dataTipTextRow('Z', app.pos(:,3), '%.3f') ...
                ];
            end

            axis(app.PlotAx, 'equal');
            grid(app.PlotAx, 'on');
            view(app.PlotAx, 3);
            title(app.PlotAx, 'Equilibrium structure');
            app.switchToView('3D');
        end

        %% ---------------- Callbacks: animation ----------------
        function PlayButtonPushed(app)
            if ~app.isBuilt
                app.showAlert('Click "Build & Solve" first.', 'Not built yet');
                return;
            end

            modeIdx = app.getSelectedModeIndex();
            Natoms = size(app.pos, 1);

            try
                app.modeShape = extractModeShape(app.V(:, modeIdx), app.masses, Natoms);
            catch ME
                app.showAlert(ME.message, 'Mode Extraction Error');
                return;
            end

            app.currentOmega = app.omega(modeIdx);

            if app.currentOmega < 1e-4
                app.setStatus('Selected mode has zero frequency (rigid body motion).', [0.8 0.4 0.0]);
            end

            app.setupAnimationHandles();
            app.animTime = 0;

            period = 1 / app.SpeedSlider.Value;
            app.animTimer = timer('ExecutionMode', 'fixedRate', 'Period', period, ...
                'TimerFcn', @(~, ~) app.animationTick(), ...
                'ErrorFcn', @(~, ~) app.stopTimerIfRunning());

            start(app.animTimer);

            app.PlayButton.Enable = 'off';
            app.PauseButton.Enable = 'on';
            app.StopButton.Enable = 'on';
        end

        function setupAnimationHandles(app)
            app.ensurePlotFigure();
            cla(app.PlotAx);
            hold(app.PlotAx, 'on');

            Nbonds = size(app.bonds, 1);
            app.bondHandles = gobjects(Nbonds, 1);
            for kk = 1:Nbonds
                i = app.bonds(kk, 1);
                j = app.bonds(kk, 2);
                app.bondHandles(kk) = plot3(app.PlotAx, ...
                    [app.pos(i,1) app.pos(j,1)], [app.pos(i,2) app.pos(j,2)], [app.pos(i,3) app.pos(j,3)], ...
                    'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
            end

            % Map elements to colors and names
            [atomColors, elementNames] = getElementDetails(app.masses);

            app.atomHandle = scatter3(app.PlotAx, app.pos(:,1), app.pos(:,2), app.pos(:,3), ...
                160, atomColors, 'filled', 'MarkerEdgeColor', 'k');

            % Add hover / click data tips safely (R2019b+)
            if isprop(app.atomHandle, 'DataTipTemplate') && exist('dataTipTextRow', 'file')
                dt = app.atomHandle.DataTipTemplate;
                dt.DataTipRows = [ ...
                    dataTipTextRow('Element', elementNames); ...
                    dataTipTextRow('Mass', app.masses, '%.3f amu'); ...
                    dataTipTextRow('X', app.pos(:,1), '%.3f'); ...
                    dataTipTextRow('Y', app.pos(:,2), '%.3f'); ...
                    dataTipTextRow('Z', app.pos(:,3), '%.3f') ...
                ];
            end

            maxSwing = app.AmplitudeSlider.Value * max(sqrt(sum(app.modeShape.^2, 2)));
            margin = maxSwing * 1.3 + 0.15;
            lo = min(app.pos, [], 1) - margin;
            hi = max(app.pos, [], 1) + margin;
            xlim(app.PlotAx, [lo(1) hi(1)]);
            ylim(app.PlotAx, [lo(2) hi(2)]);
            zlim(app.PlotAx, [lo(3) hi(3)]);
            axis(app.PlotAx, 'equal');
            grid(app.PlotAx, 'on');
            view(app.PlotAx, 3);
            title(app.PlotAx, sprintf('Animating mode %d (\\omega = %.4f)', ...
                app.getSelectedModeIndex(), app.currentOmega));
            app.switchToView('3D');
        end

        function animationTick(app)
            if isempty(app.PlotFig) || ~isvalid(app.PlotFig)
                app.stopTimerIfRunning();
                return;
            end
            app.animTime = app.animTime + app.animTimer.Period;
            amp = app.AmplitudeSlider.Value;
            displaced = app.pos + amp * sin(app.currentOmega * app.animTime) * app.modeShape;

            set(app.atomHandle, 'XData', displaced(:,1), 'YData', displaced(:,2), 'ZData', displaced(:,3));
            for kk = 1:size(app.bonds, 1)
                i = app.bonds(kk, 1);
                j = app.bonds(kk, 2);
                set(app.bondHandles(kk), ...
                    'XData', [displaced(i,1) displaced(j,1)], ...
                    'YData', [displaced(i,2) displaced(j,2)], ...
                    'ZData', [displaced(i,3) displaced(j,3)]);
            end
            drawnow limitrate;
        end

        function PauseButtonPushed(app)
            app.stopTimerIfRunning();
            app.PlayButton.Enable = 'on';
            app.PauseButton.Enable = 'off';
        end

        function StopButtonPushed(app)
            app.stopTimerIfRunning();
            app.plotStaticLattice();
            app.PlayButton.Enable = 'on';
            app.PauseButton.Enable = 'off';
            app.StopButton.Enable = 'off';
        end

        function SpeedSliderValueChanged(app)
            if ~isempty(app.animTimer) && isvalid(app.animTimer) && strcmp(app.animTimer.Running, 'on')
                stop(app.animTimer);
                app.animTimer.Period = 1 / app.SpeedSlider.Value;
                start(app.animTimer);
            end
        end

        %% ---------------- Small helpers ----------------
        function setStatus(app, msg, color)
            app.StatusLabel.Text = msg;
            if nargin > 2
                app.StatusLabel.FontColor = color;
            end
            drawnow;
        end

        function showAlert(app, msg, titleStr)
            if exist('uialert', 'file') || exist('uialert', 'builtin')
                uialert(app.UIFigure, msg, titleStr);
            else
                errordlg(msg, titleStr);
            end
        end

        function stopTimerIfRunning(app)
            if ~isempty(app.animTimer) && isvalid(app.animTimer)
                stop(app.animTimer);
                delete(app.animTimer);
            end
        end

        function onClose(app)
            app.stopTimerIfRunning();
            if ~isempty(app.PlotFig) && isvalid(app.PlotFig)
                delete(app.PlotFig);
            end
            delete(app.UIFigure);
        end
    end
end