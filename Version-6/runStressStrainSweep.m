%% runStressStrainSweep.m
% Standalone script: sweeps applied engineering shear strain across a
% mass-spring crystal lattice, tracks average stress and elastic energy
% at every step, extracts the shear modulus G from the linear-elastic
% regime, and marks the yield point predicted by the Von Mises / Tresca
% criteria.
%
% This builds directly on the vibrational engine you already have
% (generateLatticeGeneral.m, buildDynamicalMatrix.m, jacobiEigenSolver.m)
% plus the new static-load layer (computeShearForces.m,
% computeAtomicStress.m, checkPlasticity.m).
%
% Study Guide reference: Sec. 2-4 (deformation/force/stress/yield) and
% Sec. 6 (elastic modulus extraction, stress-strain curve).

clear; clc; close all;

%% ---- 1. Build the baseline lattice ------------------------------------
% A simple cubic lattice is used here so this script is fully
% self-contained and runnable with no external .cif file. Swap this
% block for parseCIFFile.m + generateLatticeGeneral.m's CIF-based
% calling form if you want to run the sweep on one of your real
% crystal structures instead.
a0 = 1.0;                                  % lattice constant (arb. units)
latticeVectors = a0 * eye(3);              % simple cubic primitive cell
basisFrac      = [0 0 0];                  % one atom per primitive cell
basisMasses    = 1.0;                      % arbitrary atomic mass
N         = 4;                             % N x N x N supercell
numShells = 1;                             % nearest-neighbor bonds only

[pos, bonds, masses] = generateLatticeGeneral(latticeVectors, basisFrac, basisMasses, N, numShells);
Natoms = size(pos, 1);

springConstant = 1.0;                      % k (arb. units)
K_matrix = buildDynamicalMatrix(pos, bonds, springConstant, masses);

% Atomic volume Omega_i (Study Guide Sec. 3.1): unit-cell volume divided
% by the number of basis atoms per primitive cell. Every atom is
% equivalent here (monatomic simple cubic), so a single scalar suffices;
% computeAtomicStress.m also accepts a per-atom [Natoms x 1] vector for
% multi-species/multi-site lattices.
cellVolume   = abs(det(latticeVectors));
atomsPerCell = size(basisFrac, 1);
atomicVolume = cellVolume / atomsPerCell;

%% ---- 2. Yield strength inputs ------------------------------------------
% Illustrative round-number thresholds in the same (arbitrary) stress
% units as springConstant/a0^2. Replace with literature values (e.g.
% Callister, Ch. 6) if you calibrate this lattice to a real material.
sigmaYield = 0.15;      % Von Mises yield strength
tauYield   = 0.10;      % Tresca (max shear) yield strength

%% ---- 3. Sweep the applied engineering shear strain gamma ---------------
gammaMax = 0.30;
nSteps   = 60;
gammaVec = linspace(0, gammaMax, nSteps)';

avgVonMises   = zeros(nSteps, 1);
avgSigmaXY    = zeros(nSteps, 1);   % macroscopic shear stress tau_xy
elasticEnergy = zeros(nSteps, 1);
anyYieldVM    = false(nSteps, 1);
anyYieldTr    = false(nSteps, 1);

for s = 1:nSteps
    gamma = gammaVec(s);

    % --- Step 1: affine shear + bond-by-bond Hooke's law (Sec. 2) ---
    [~, dispField, bondVecs_sh, ~, bondForceVec, ~] = ...
        computeShearForces(pos, bonds, springConstant, gamma);

    % --- Step 2: per-atom virial stress + principal stresses (Sec. 3) ---
    [stressTensor, principalStresses, vonMises] = ...
        computeAtomicStress(bonds, bondVecs_sh, bondForceVec, Natoms, atomicVolume);

    % --- Step 3: yield criteria (Sec. 4) ---
    [vmYield, trYield] = checkPlasticity(principalStresses, vonMises, sigmaYield, tauYield);
    anyYieldVM(s) = any(vmYield);
    anyYieldTr(s) = any(trYield);

    % --- Step 4: macroscopic response for the stress-strain curve ---
    avgVonMises(s) = mean(vonMises);
    avgSigmaXY(s)  = mean(squeeze(stressTensor(1, 2, :)));   % mean sigma_xy

    % Elastic energy E = 1/2 * u' * K * u (Sec. 6.2), reusing the SAME
    % stiffness matrix K already built for the vibrational problem.
    % dispField is [Natoms x 3]; K's DOF ordering is
    % [atom1_x, atom1_y, atom1_z, atom2_x, ...] (see
    % buildDynamicalMatrix.m: idx_i = 3*(i-1)+(1:3)), so transposing
    % dispField before flattening puts u into that exact same order.
    u = reshape(dispField', 3 * Natoms, 1);
    elasticEnergy(s) = 0.5 * (u' * K_matrix * u);
end

%% ---- 4. Extract the shear modulus G from the linear-elastic regime -----
% tau_xy = G*gamma is the continuum definition of the shear modulus
% (Study Guide Sec. 2.4 / 6.1); fit it on the pre-yield points only,
% using the FIRST index (if any) where either criterion trips.
firstYieldIdx = find(anyYieldVM | anyYieldTr, 1);
if isempty(firstYieldIdx)
    linearIdx = 1:nSteps;                       % never yielded in this sweep
else
    linearIdx = 1:max(firstYieldIdx - 1, 2);
end

pFit  = polyfit(gammaVec(linearIdx), avgSigmaXY(linearIdx), 1);
G_fit = pFit(1);
fprintf('Fitted shear modulus G = %.6g (from tau_xy = G*gamma, %d-point linear fit)\n', ...
    G_fit, numel(linearIdx));

% Frenkel ideal shear strength (Study Guide Sec. 5), for reference only -
% real crystals yield 100-1000x below this due to dislocation glide.
tau_ideal = G_fit / (2 * pi);
fprintf('Frenkel ideal shear strength estimate: tau_ideal ~= G/(2*pi) = %.6g\n', tau_ideal);

%% ---- 5. Plot: stress-strain curve + elastic energy ----------------------
% Formatting: thick lines, large fonts, high-contrast colors, suitable
% for direct export into a LaTeX document (e.g. via exportgraphics).
lw = 2; fs = 13;

figure('Color', 'w', 'Position', [100 100 950 750]);

% --- Panel 1: stress-strain curve ---
subplot(2,1,1); hold on; box on; grid on;
plot(gammaVec, avgSigmaXY, '-o', 'LineWidth', lw, 'Color', [0.00 0.45 0.74], ...
    'MarkerFaceColor', [0.00 0.45 0.74], 'MarkerSize', 5, ...
    'DisplayName', '\tau_{xy} (mean shear stress)');
plot(gammaVec, avgVonMises, '-s', 'LineWidth', lw, 'Color', [0.85 0.33 0.10], ...
    'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize', 5, ...
    'DisplayName', '\sigma_{vm} (mean von Mises)');
plot(gammaVec(linearIdx), polyval(pFit, gammaVec(linearIdx)), '--', ...
    'LineWidth', lw, 'Color', [0.2 0.2 0.2], ...
    'DisplayName', sprintf('linear fit, G = %.4g', G_fit));
yline(sigmaYield, ':', '\sigma_{yield}', 'LineWidth', lw, 'Color', [0.6 0 0], ...
    'FontSize', fs, 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
if ~isempty(firstYieldIdx)
    xline(gammaVec(firstYieldIdx), '-', sprintf('yield at \\gamma = %.3f', gammaVec(firstYieldIdx)), ...
        'LineWidth', lw, 'Color', [0.6 0 0], 'FontSize', fs, ...
        'LabelOrientation', 'horizontal', 'HandleVisibility', 'off');
    plot(gammaVec(firstYieldIdx), avgVonMises(firstYieldIdx), 'p', ...
        'MarkerSize', 16, 'MarkerFaceColor', [1 0.6 0], 'MarkerEdgeColor', 'k', ...
        'HandleVisibility', 'off');
end
xlabel('Engineering shear strain \gamma', 'FontSize', fs);
ylabel('Stress', 'FontSize', fs);
title('Stress-Strain Curve and Yield Point', 'FontSize', fs + 2, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', fs - 2);
set(gca, 'FontSize', fs, 'LineWidth', 1.2);

% --- Panel 2: elastic energy curve ---
subplot(2,1,2); hold on; box on; grid on;
plot(gammaVec, elasticEnergy, '-^', 'LineWidth', lw, 'Color', [0.47 0.67 0.19], ...
    'MarkerFaceColor', [0.47 0.67 0.19], 'MarkerSize', 5);
if ~isempty(firstYieldIdx)
    xline(gammaVec(firstYieldIdx), '-', 'LineWidth', lw, 'Color', [0.6 0 0]);
    plot(gammaVec(firstYieldIdx), elasticEnergy(firstYieldIdx), 'p', ...
        'MarkerSize', 16, 'MarkerFaceColor', [1 0.6 0], 'MarkerEdgeColor', 'k');
end
xlabel('Engineering shear strain \gamma', 'FontSize', fs);
ylabel('Elastic energy  E = 1/2 u^{T}Ku', 'FontSize', fs);
title('Elastic Strain Energy', 'FontSize', fs + 2, 'FontWeight', 'bold');
set(gca, 'FontSize', fs, 'LineWidth', 1.2);

sgtitle('Static Shear Loading of Crystal Lattice: Elasticity to Yield', ...
    'FontSize', fs + 4, 'FontWeight', 'bold');

% NOTE (Study Guide Sec. 6.3): because this spring lattice is purely
% harmonic, Hooke's law has no built-in upper bound - the curve above
% will keep looking smooth/roughly-linear PAST the marked yield point
% too. The star marker means "the yield criterion has been exceeded",
% NOT "the simulation itself shows plastic flow" - a defect-free spring
% lattice with no dislocations, vacancies, or grain boundaries is
% structurally incapable of representing dislocation-mediated plastic
% deformation (Study Guide Sec. 5). State this explicitly as a model
% limitation in your report.
