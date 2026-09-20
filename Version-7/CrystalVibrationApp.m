classdef CrystalVibrationApp < handle
    %CRYSTALVIBRATIONAPP Interactive crystal vibration viewer with element
    %   colors, data tips, and shear analysis - all in a SINGLE window.
    %
    %   Layout: one uifigure holding a uigridlayout with a fixed-width
    %   control column on the left and a tab group on the right
    %   ("Vibration" = 3D lattice + frequency spectrum side by side,
    %   "Shear Analysis" = pristine/sheared lattices and their spectra).
    %   Every control is placed in a grid cell rather than at absolute
    %   pixel coordinates, so nothing can overlap at any window size.

    properties (Access = private)
        % --- Single application window ---
        UIFigure
        MainGrid
        ControlPanel
        TabGroup
        VibrationTab
        ShearTab

        % --- Model / animation controls (left column) ---
        LoadFileButton
        StatusLabel
        NEditField
        kEditField
        cutoffDropDown
        BuildButton
        ModeDropDown
        AmplitudeSlider
        SpeedSlider
        PlayPauseButton
        StopButton

        % --- Vibration tab ---
        PlotAx
        SpecAx
        SpecHighlightHandle

        % --- Shear Analysis tab ---
        ShearSlider
        RelaxButton
        RelaxStatusLabel
        UIAxes3D_NoForce
        UIAxes3D_WithForce
        UIAxesMode_NoForce
        UIAxesMode_WithForce
        shearRendered = false   % shear tab is rendered lazily - it re-solves the eigenproblem

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
            if isvalid(app.UIFigure)
                delete(app.UIFigure);
            end
        end
    end

    methods (Access = private)

        %% ---------------- UI construction ----------------
        function createComponents(app)
            app.UIFigure = uifigure('Name', 'Crystal Vibration Explorer', ...
                'Position', [80 80 1360 820], 'Color', 'w');
            app.UIFigure.CloseRequestFcn = @(~, ~) app.onClose();

            app.MainGrid = uigridlayout(app.UIFigure, [1 2]);
            app.MainGrid.ColumnWidth = {330, '1x'};
            app.MainGrid.RowHeight = {'1x'};
            app.MainGrid.Padding = [10 10 10 10];
            app.MainGrid.ColumnSpacing = 10;

            app.buildControlPanel();

            app.TabGroup = uitabgroup(app.MainGrid);
            app.TabGroup.Layout.Row = 1;
            app.TabGroup.Layout.Column = 2;
            app.TabGroup.SelectionChangedFcn = @(~, ~) app.TabSelectionChanged();

            app.VibrationTab = uitab(app.TabGroup, 'Title', 'Vibration');
            app.ShearTab = uitab(app.TabGroup, 'Title', 'Shear Analysis');

            app.buildVibrationTab();
            app.buildShearTab();
        end

        function buildControlPanel(app)
            app.ControlPanel = uipanel(app.MainGrid, ...
                'Title', 'Model & Animation Controls', ...
                'FontWeight', 'bold', 'BackgroundColor', 'w');
            app.ControlPanel.Layout.Row = 1;
            app.ControlPanel.Layout.Column = 1;

            g = uigridlayout(app.ControlPanel, [17 1]);
            g.ColumnWidth = {'1x'};
            g.RowHeight = {32, 34, 20, 25, 20, 25, 20, 25, 38, 20, 25, 20, 55, 20, 55, 34, '1x'};
            g.Padding = [10 10 10 10];
            g.RowSpacing = 6;

            app.LoadFileButton = uibutton(g, 'push', 'Text', 'Load Crystal File...', ...
                'ButtonPushedFcn', @(~, ~) app.LoadFileButtonPushed());
            app.LoadFileButton.Layout.Row = 1;

            app.StatusLabel = uilabel(g, 'Text', 'No file loaded (using default bcc structure)', ...
                'FontColor', [0.4 0.4 0.4], 'WordWrap', 'on');
            app.StatusLabel.Layout.Row = 2;

            lbl = uilabel(g, 'Text', 'Supercell size N (NxNxN cells):');
            lbl.Layout.Row = 3;
            app.NEditField = uieditfield(g, 'numeric', 'Value', 2, 'Limits', [1 4]);
            app.NEditField.Layout.Row = 4;

            lbl = uilabel(g, 'Text', 'Spring constant k:');
            lbl.Layout.Row = 5;
            app.kEditField = uieditfield(g, 'numeric', 'Value', 1.0, 'Limits', [0.001 Inf]);
            app.kEditField.Layout.Row = 6;

            lbl = uilabel(g, 'Text', 'Bond network:');
            lbl.Layout.Row = 7;
            app.cutoffDropDown = uidropdown(g, ...
                'Items', {'Nearest neighbors only', '1st + 2nd neighbors (recommended, rigid)'}, ...
                'Value', '1st + 2nd neighbors (recommended, rigid)');
            app.cutoffDropDown.Layout.Row = 8;

            app.BuildButton = uibutton(g, 'push', 'Text', 'Build && Solve', ...
                'FontWeight', 'bold', 'BackgroundColor', [0.20 0.45 0.80], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.BuildButtonPushed());
            app.BuildButton.Layout.Row = 9;

            lbl = uilabel(g, 'Text', 'Vibrational mode:');
            lbl.Layout.Row = 10;
            app.ModeDropDown = uidropdown(g, 'Items', {'(build the model first)'}, ...
                'Enable', 'off', 'ValueChangedFcn', @(~, ~) app.ModeDropDownValueChanged());
            app.ModeDropDown.Layout.Row = 11;

            lbl = uilabel(g, 'Text', 'Amplitude:');
            lbl.Layout.Row = 12;
            app.AmplitudeSlider = uislider(g, 'Limits', [0 1], 'Value', 0.25);
            app.AmplitudeSlider.Layout.Row = 13;

            lbl = uilabel(g, 'Text', 'Speed (fps):');
            lbl.Layout.Row = 14;
            app.SpeedSlider = uislider(g, 'Limits', [5 60], 'Value', 30, ...
                'ValueChangedFcn', @(~, ~) app.SpeedSliderValueChanged());
            app.SpeedSlider.Layout.Row = 15;

            % Play and Pause were mutually exclusive (one was always
            % disabled), so they are a single toggle here; Stop stays
            % separate because it also resets to the rest positions.
            btnGrid = uigridlayout(g, [1 2]);
            btnGrid.Layout.Row = 16;
            btnGrid.ColumnWidth = {'1x', '1x'};
            btnGrid.Padding = [0 0 0 0];
            btnGrid.ColumnSpacing = 6;

            app.PlayPauseButton = uibutton(btnGrid, 'push', 'Text', 'Play', ...
                'Enable', 'off', 'BackgroundColor', [0.30 0.70 0.35], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.PlayPauseButtonPushed());
            app.PlayPauseButton.Layout.Column = 1;

            app.StopButton = uibutton(btnGrid, 'push', 'Text', 'Stop', 'Enable', 'off', ...
                'ButtonPushedFcn', @(~, ~) app.StopButtonPushed());
            app.StopButton.Layout.Column = 2;
        end

        function buildVibrationTab(app)
            g = uigridlayout(app.VibrationTab, [1 2]);
            g.ColumnWidth = {'1.25x', '1x'};
            g.RowHeight = {'1x'};
            g.Padding = [8 8 8 8];
            g.ColumnSpacing = 8;

            app.PlotAx = uiaxes(g);
            app.PlotAx.Layout.Row = 1;
            app.PlotAx.Layout.Column = 1;
            title(app.PlotAx, 'Click "Build && Solve" to begin');
            app.PlotAx.DataAspectRatio = [1 1 1];
            grid(app.PlotAx, 'on');
            view(app.PlotAx, 3);

            app.SpecAx = uiaxes(g);
            app.SpecAx.Layout.Row = 1;
            app.SpecAx.Layout.Column = 2;
            title(app.SpecAx, 'Frequency spectrum');
            xlabel(app.SpecAx, 'Mode index');
            ylabel(app.SpecAx, '\omega (angular frequency)');
            grid(app.SpecAx, 'on');
        end

        function buildShearTab(app)
            g = uigridlayout(app.ShearTab, [3 2]);
            g.ColumnWidth = {'1x', '1x'};
            g.RowHeight = {90, '1x', '1x'};
            g.Padding = [8 8 8 8];
            g.RowSpacing = 8;
            g.ColumnSpacing = 8;

            ctrl = uigridlayout(g, [2 3]);
            ctrl.Layout.Row = 1;
            ctrl.Layout.Column = [1 2];
            ctrl.ColumnWidth = {'1x', 220, '1.1x'};
            ctrl.RowHeight = {20, '1x'};
            ctrl.Padding = [0 0 0 0];
            ctrl.ColumnSpacing = 12;

            lbl = uilabel(ctrl, 'Text', 'Shear force magnitude:');
            lbl.Layout.Row = 1;
            lbl.Layout.Column = 1;

            app.ShearSlider = uislider(ctrl, 'Limits', [0 2], 'Value', 0, ...
                'ValueChangedFcn', @(src, ~) app.updateShearAnalysis(src.Value));
            app.ShearSlider.Layout.Row = 2;
            app.ShearSlider.Layout.Column = 1;

            % Relaxation is a one-shot linear solve (relaxShearLattice.m),
            % too slow to run on every slider tick for a larger supercell,
            % so it stays an on-demand button rather than a slider callback.
            app.RelaxButton = uibutton(ctrl, 'push', 'Text', 'Relax (Gauss-Jordan, Exp.5)', ...
                'FontWeight', 'bold', 'BackgroundColor', [0.20 0.45 0.80], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~, ~) app.RelaxButtonPushed());
            app.RelaxButton.Layout.Row = 2;
            app.RelaxButton.Layout.Column = 2;

            app.RelaxStatusLabel = uilabel(ctrl, ...
                'Text', 'Build & Solve first, then drag the slider. Relax solves for true equilibrium of the interior atoms.', ...
                'FontColor', [0.4 0.4 0.4], 'WordWrap', 'on');
            app.RelaxStatusLabel.Layout.Row = [1 2];
            app.RelaxStatusLabel.Layout.Column = 3;

            app.UIAxes3D_NoForce = uiaxes(g);
            app.UIAxes3D_NoForce.Layout.Row = 2;
            app.UIAxes3D_NoForce.Layout.Column = 1;
            title(app.UIAxes3D_NoForce, 'Pristine Lattice (No Force)');

            app.UIAxes3D_WithForce = uiaxes(g);
            app.UIAxes3D_WithForce.Layout.Row = 2;
            app.UIAxes3D_WithForce.Layout.Column = 2;
            title(app.UIAxes3D_WithForce, 'Shear-Deformed Lattice (Force Distributed)');

            app.UIAxesMode_NoForce = uiaxes(g);
            app.UIAxesMode_NoForce.Layout.Row = 3;
            app.UIAxesMode_NoForce.Layout.Column = 1;
            title(app.UIAxesMode_NoForce, 'Normal Modes - Pristine');
            xlabel(app.UIAxesMode_NoForce, 'Mode Index');
            ylabel(app.UIAxesMode_NoForce, '\omega (angular frequency)');

            app.UIAxesMode_WithForce = uiaxes(g);
            app.UIAxesMode_WithForce.Layout.Row = 3;
            app.UIAxesMode_WithForce.Layout.Column = 2;
            title(app.UIAxesMode_WithForce, 'Normal Modes - Sheared');
            xlabel(app.UIAxesMode_WithForce, 'Mode Index');
            ylabel(app.UIAxesMode_WithForce, '\omega (angular frequency)');
        end

        function TabSelectionChanged(app)
            % The shear tab re-solves the eigenproblem twice, so only
            % render it when it is actually brought to the front.
            if app.TabGroup.SelectedTab == app.ShearTab && app.isBuilt && ~app.shearRendered
                app.updateShearAnalysis(app.ShearSlider.Value);
            end
        end

        %% ---------------- Shear analysis ----------------
        function updateShearAnalysis(app, shearVal)
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
            title(app.UIAxes3D_NoForce, 'Pristine Lattice (No Force)');

            cla(app.UIAxes3D_WithForce);
            plotLatticeWithGradient(app.UIAxes3D_WithForce, coords_sh, elementNames, force_dist);
            title(app.UIAxes3D_WithForce, 'Shear-Deformed Lattice (Force Distributed)');

            cla(app.UIAxesMode_NoForce);
            plot(app.UIAxesMode_NoForce, freqs_base, '-o', 'LineWidth', 1.5, 'Color', 'b');

            cla(app.UIAxesMode_WithForce);
            plot(app.UIAxesMode_WithForce, freqs_sh, '-o', 'LineWidth', 1.5, 'Color', 'r');
            title(app.UIAxesMode_WithForce, 'Normal Modes - Sheared');

            app.shearRendered = true;
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
            app.setPlayPauseState('play');
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
                app.PlayPauseButton.Enable = 'on';
                app.StopButton.Enable = 'off';

                % The shear tab is now stale; re-render it if it is the
                % visible tab, otherwise leave it for TabSelectionChanged.
                app.shearRendered = false;
                if app.TabGroup.SelectedTab == app.ShearTab
                    app.updateShearAnalysis(app.ShearSlider.Value);
                end
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
            title(app.SpecAx, 'Frequency spectrum');
            hold(app.SpecAx, 'off');
            app.SpecHighlightHandle = []; % axis was cleared, so any old marker handle is stale
        end

        %% ---------------- Static 3D lattice plot (with element colors + data tips) ----------------
        function plotStaticLattice(app)
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

            % Data tips: element, mass, and atom index. Data tips are
            % available on UIAxes by default (click an atom), so no
            % datacursormode call is needed - and leaving cursor mode off
            % keeps click-drag free for rotating the 3D view.
            app.atomHandle.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Element', elementNames);
            app.atomHandle.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Mass', app.masses);
            app.atomHandle.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Atom #', (1:numel(app.masses))');

            ax.DataAspectRatio = [1 1 1];
            grid(ax, 'on');
            view(ax, 3);
            title(ax, sprintf('%d atoms, %d bonds (click an atom for details)', size(app.pos, 1), nBonds));
            hold(ax, 'off');
        end

        %% ---------------- Mode spectrum highlight ----------------
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
        function PlayPauseButtonPushed(app)
            if ~app.isBuilt
                return;
            end

            if ~isempty(app.animTimer) && isvalid(app.animTimer) && strcmp(app.animTimer.Running, 'on')
                stop(app.animTimer);
                app.setPlayPauseState('play');
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

            app.setPlayPauseState('pause');
            app.StopButton.Enable = 'on';
        end

        function setPlayPauseState(app, state)
            % 'play'  -> button offers Play (animation is stopped/paused)
            % 'pause' -> button offers Pause (animation is running)
            if strcmp(state, 'pause')
                app.PlayPauseButton.Text = 'Pause';
                app.PlayPauseButton.BackgroundColor = [0.85 0.65 0.15];
            else
                app.PlayPauseButton.Text = 'Play';
                app.PlayPauseButton.BackgroundColor = [0.30 0.70 0.35];
            end
        end

        function StopButtonPushed(app)
            app.stopTimerIfRunning();
            app.animTime = 0;
            if app.isBuilt
                app.plotStaticLattice(); % snap back to rest positions
            end
            app.setPlayPauseState('play');
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

% applyShearForce, plotLattice, plotLatticeWithGradient, and
% relaxShearLattice are provided as separate files on the MATLAB path.
% They are intentionally NOT redefined here as local functions: a local
% function of the same name in this class file would shadow the real
% path versions for every call made from within this class, silently
% substituting the wrong implementation. Make sure those files sit next
% to this one (or anywhere else on the path) alongside
% buildDynamicalMatrix.m, jacobiEigenSolver.m, parseCIFFile.m,
% parseCrystalFile.m, gaussJordanSolve.m, and getBuiltinCrystalDef.m,
% none of which are defined in this file either.
