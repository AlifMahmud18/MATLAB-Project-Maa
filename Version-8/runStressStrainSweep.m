%% runStressStrainSweep.m
% Standalone script: sweeps applied engineering shear strain across a
% mass-spring crystal lattice, tracks average stress and elastic energy
% at every step, extracts the shear modulus G from the linear-elastic
% regime, and pinpoints the yield strain predicted by the Von Mises /
% Tresca criteria.
%
% This builds directly on the vibrational engine you already have
% (generateLatticeGeneral.m, buildDynamicalMatrix.m, jacobiEigenSolver.m)
% plus the static-load layer (computeShearForces.m, computeAtomicStress.m,
% checkPlasticity.m), and now routes every post-processing step through
% the EEE 212 numerical-methods toolkit instead of MATLAB's built-ins:
%
%   Exp. 4  polyLeastSquaresFit.m   - shear-modulus fit (was polyfit)
%   Exp. 5  gaussJordanSolve.m      - linear solve inside the fit above
%   Exp. 7  compositeTrapz.m        - strain-energy from the stress curve
%   Exp. 8  bisectionRoot.m /       - exact yield strain (was just the
%           newtonRaphsonRoot.m       nearest coarse sample)
%   Exp. 3  newtonDividedDiffInterp.m - interpolated cross-check at an
%                                       off-grid strain value
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
numShells = 2;                             % 1st + 2nd neighbor bonds - nearest-neighbor-only
                                            % leaves this lattice severely under-constrained
                                            % (Maxwell rigidity: rank(K) came out 144/192, i.e.
                                            % 48 floppy zero-energy modes instead of the 3 expected
                                            % rigid-body translations), which was corrupting the
                                            % shear response below. This matches the GUI's own
                                            % "1st + 2nd neighbors (recommended, rigid)" default.

[pos, bonds, masses] = generateLatticeGeneral(latticeVectors, basisFrac, basisMasses, N, numShells);
Natoms = size(pos, 1);

springConstant = 1.0;                      % k (arb. units)

% Atomic volume Omega_i (Study Guide Sec. 3.1): unit-cell volume divided
% by the number of basis atoms per primitive cell. Every atom is
% equivalent here (monatomic simple cubic), so a single scalar suffices;
% computeAtomicStress.m also accepts a per-atom [Natoms x 1] vector for
% multi-species/multi-site lattices.
cellVolume   = abs(det(latticeVectors));
atomsPerCell = size(basisFrac, 1);
atomicVolume = cellVolume / atomsPerCell;
totalVolume  = cellVolume * N^3;           % whole supercell, for the Exp.7 energy check

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
h_gamma  = gammaVec(2) - gammaVec(1);       % uniform step, needed by compositeTrapz.m

avgVonMises   = zeros(nSteps, 1);
avgSigmaXY    = zeros(nSteps, 1);   % macroscopic shear stress tau_xy
elasticEnergy = zeros(nSteps, 1);
anyYieldVM    = false(nSteps, 1);
anyYieldTr    = false(nSteps, 1);

for s = 1:nSteps
    gamma = gammaVec(s);

    % --- Step 1: affine shear + bond-by-bond Hooke's law (Sec. 2) ---
    [~, ~, bondVecs_sh, bondForceMag, bondForceVec, ~] = ...
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

    % Elastic (bond) energy, E = sum_bonds 1/2*k*(d_new-d0)^2, reusing the
    % SAME bondForceMag = k*(d_new-d0) already computed above (kept as
    % the analytic cross-check for the Exp.7 numerical-integration result
    % in Section 4b below).
    %
    % NOTE: this is NOT computed as 0.5*u'*K*u with u = the raw affine
    % displacement field. generateLatticeGeneral.m builds bonds under
    % periodic minimum-image wrapping, so a good fraction of bonds
    % connect an atom to a periodic IMAGE of its neighbor (37.5% of
    % bonds for this exact lattice/cutoff). K's own matrix entries
    % correctly encode each bond's true (wrapped) geometry, but a global
    % "u" array built from the raw, UNWRAPPED affine formula u=gamma*S*r
    % does not - for a wrapped bond, u(j)-u(i) computed that way uses the
    % wrong (unwrapped, often huge) separation instead of the true short
    % bond vector, silently inflating 0.5*u'*K*u (verified: this
    % overstated the energy by ~3x-8x depending on supercell size, while
    % never converging back to 1x even as gamma -> 0). bondForceMag was
    % already computed correctly above (Section 1 uses the wraparound-
    % safe bond vectors from computeShearForces.m), so reusing it here
    % sidesteps the issue entirely instead of re-deriving it incorrectly.
    elasticEnergy(s) = 0.5 * sum(bondForceMag.^2) / springConstant;
end

%% ---- 4a. Shear modulus via least-squares regression (Exp. 4 + Exp. 5) --
% tau_xy = G*gamma is the continuum definition of the shear modulus
% (Study Guide Sec. 2.4 / 6.1); fit it on the pre-yield points only,
% using the FIRST coarse sample index (if any) where either criterion
% trips - this only needs to bracket the linear regime, the exact yield
% strain is pinned down precisely in Section 4c below.
firstYieldIdx = find(anyYieldVM | anyYieldTr, 1);
if isempty(firstYieldIdx)
    linearIdx = 1:nSteps;                       % never yielded in this sweep
else
    linearIdx = 1:max(firstYieldIdx - 1, 2);
end

% polyLeastSquaresFit.m builds the same normal equations as the labsheet
% (Exp. 4) and solves them with gaussJordanSolve.m (Exp. 5) instead of
% MATLAB's polyfit/backslash. Coefficients come back ASCENDING: [a0 a1].
cFit  = polyLeastSquaresFit(gammaVec(linearIdx), avgSigmaXY(linearIdx), 1);
G_fit = cFit(2);
fprintf('Fitted shear modulus G = %.6g (from tau_xy = G*gamma, %d-point linear fit, Exp.4/Exp.5)\n', ...
    G_fit, numel(linearIdx));

% Frenkel ideal shear strength (Study Guide Sec. 5), for reference only -
% real crystals yield 100-1000x below this due to dislocation glide.
tau_ideal = G_fit / (2 * pi);
fprintf('Frenkel ideal shear strength estimate: tau_ideal ~= G/(2*pi) = %.6g\n', tau_ideal);

%% ---- 4b. Strain energy via numerical integration (Exp. 7) --------------
% Classical result: for a linear stress-strain relation the strain
% energy DENSITY is the area under the stress-strain curve,
% w(gamma) = integral_0^gamma tau_xy(g) dg, and E = w * (total volume).
% compositeTrapz.m implements the labsheet's composite trapezoidal rule
% exactly; comparing its result to the exact bond-by-bond spring energy
% computed above (Section 3) is a direct numerical-methods validation of
% the physics model - the two should agree to O(h_gamma^2).
[~, energyDensityCum] = compositeTrapz(avgSigmaXY, h_gamma);
elasticEnergy_numeric = energyDensityCum * totalVolume;

energyCheckErr = max(abs(elasticEnergy_numeric - elasticEnergy));
energyCheckRel = energyCheckErr / max(abs(elasticEnergy(end)), eps);
fprintf(['Exp.7 cross-check: trapezoidal-integrated energy vs exact bond-sum energy -> ' ...
    'max abs diff = %.3e (relative %.3e)\n'], energyCheckErr, energyCheckRel);

%% ---- 4c. Exact yield strain via root-finding (Exp. 8) -------------------
% Each criterion's "gap" function gap(gamma) = worst-atom stress -
% threshold crosses zero exactly at the true yield strain; shearYieldGap.m
% re-runs the real physics (computeShearForces + computeAtomicStress) at
% an arbitrary continuous gamma, so bisectionRoot.m and
% newtonRaphsonRoot.m (Exp. 8) can bracket/iterate to that root directly
% instead of only reporting the nearest of the nSteps coarse samples.
gammaYieldVM = NaN; gammaYieldTr = NaN;
itersBisVM = 0; itersNewtVM = 0;
itersBisTr = 0; itersNewtTr = 0;

firstYieldIdxVM = find(anyYieldVM, 1);
if ~isempty(firstYieldIdxVM) && firstYieldIdxVM > 1
    gVM = @(g) shearYieldGap(g, pos, bonds, springConstant, Natoms, atomicVolume, 'vonmises', sigmaYield);
    [gammaYieldVM, ~, itersBisVM] = bisectionRoot(gVM, gammaVec(firstYieldIdxVM-1), gammaVec(firstYieldIdxVM), 1e-8);
    [gammaYieldVM_newton, itersNewtVM] = newtonRaphsonRoot(gVM, gammaYieldVM, 1e-8);
    fprintf(['Von Mises yield strain: bisection gamma = %.8f (%d iters), ' ...
        'Newton-Raphson gamma = %.8f (%d iters), agreement = %.2e\n'], ...
        gammaYieldVM, itersBisVM, gammaYieldVM_newton, itersNewtVM, abs(gammaYieldVM - gammaYieldVM_newton));
end

firstYieldIdxTr = find(anyYieldTr, 1);
if ~isempty(firstYieldIdxTr) && firstYieldIdxTr > 1
    gTr = @(g) shearYieldGap(g, pos, bonds, springConstant, Natoms, atomicVolume, 'tresca', tauYield);
    [gammaYieldTr, ~, itersBisTr] = bisectionRoot(gTr, gammaVec(firstYieldIdxTr-1), gammaVec(firstYieldIdxTr), 1e-8);
    [gammaYieldTr_newton, itersNewtTr] = newtonRaphsonRoot(gTr, gammaYieldTr, 1e-8);
    fprintf(['Tresca yield strain:    bisection gamma = %.8f (%d iters), ' ...
        'Newton-Raphson gamma = %.8f (%d iters), agreement = %.2e\n'], ...
        gammaYieldTr, itersBisTr, gammaYieldTr_newton, itersNewtTr, abs(gammaYieldTr - gammaYieldTr_newton));
end

% Whichever criterion trips first (smaller strain) is the one that
% actually governs; NaN-safe min so a criterion that never triggers in
% this sweep doesn't win by default.
candidates = [gammaYieldVM, gammaYieldTr];
if all(isnan(candidates))
    gammaYieldRefined = [];
    yieldLabel = '';
else
    [gammaYieldRefined, which] = min(candidates);
    if which == 1
        yieldLabel = sprintf('von Mises, \\gamma = %.5f', gammaYieldRefined);
    else
        yieldLabel = sprintf('Tresca, \\gamma = %.5f', gammaYieldRefined);
    end
    fprintf('Governing (first) yield criterion: %s\n', strrep(yieldLabel, '\\', ''));
end

%% ---- 4d. Newton polynomial interpolation cross-check (Exp. 3) ----------
% Pick a strain value that does NOT land on the sweep grid, interpolate
% the tabulated (gammaVec, avgVonMises) curve there with a local Newton
% divided-difference polynomial, and compare against re-running the real
% physics at that exact strain - the same "what is f(2.5)?" idea as the
% labsheet's own Table 1 example.
if nSteps >= 6
    queryStepIdx = round(nSteps / 3);
    gammaQuery = 0.5 * (gammaVec(queryStepIdx) + gammaVec(queryStepIdx + 1));
    localIdx = max(1, queryStepIdx - 2):min(nSteps, queryStepIdx + 3);

    vmInterp = newtonDividedDiffInterp(gammaVec(localIdx), avgVonMises(localIdx), gammaQuery);

    [~, ~, bondVecs_q, ~, bondForceVec_q, ~] = computeShearForces(pos, bonds, springConstant, gammaQuery);
    [~, ~, vonMises_q] = computeAtomicStress(bonds, bondVecs_q, bondForceVec_q, Natoms, atomicVolume);
    vmTrue = mean(vonMises_q);

    fprintf(['Exp.3 interpolation check at off-grid gamma = %.5f: ' ...
        'Newton-poly estimate = %.6f, true (re-evaluated) = %.6f, error = %.2e\n'], ...
        gammaQuery, vmInterp, vmTrue, abs(vmInterp - vmTrue));
end

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
plot(gammaVec(linearIdx), cFit(1) + cFit(2) * gammaVec(linearIdx), '--', ...
    'LineWidth', lw, 'Color', [0.2 0.2 0.2], ...
    'DisplayName', sprintf('linear fit (Exp.4), G = %.4g', G_fit));
yline(sigmaYield, ':', '\sigma_{yield}', 'LineWidth', lw, 'Color', [0.6 0 0], ...
    'FontSize', fs, 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
if ~isempty(gammaYieldRefined)
    xline(gammaYieldRefined, '-', sprintf('yield (Exp.8) at %s', yieldLabel), ...
        'LineWidth', lw, 'Color', [0.6 0 0], 'FontSize', fs, ...
        'LabelOrientation', 'horizontal', 'HandleVisibility', 'off');
    plot(gammaYieldRefined, interp1(gammaVec, avgVonMises, gammaYieldRefined), 'p', ...
        'MarkerSize', 16, 'MarkerFaceColor', [1 0.6 0], 'MarkerEdgeColor', 'k', ...
        'HandleVisibility', 'off');
elseif ~isempty(firstYieldIdx)
    xline(gammaVec(firstYieldIdx), '-', sprintf('yield at \\gamma = %.3f (coarse sample)', gammaVec(firstYieldIdx)), ...
        'LineWidth', lw, 'Color', [0.6 0 0], 'FontSize', fs, ...
        'LabelOrientation', 'horizontal', 'HandleVisibility', 'off');
end
xlabel('Engineering shear strain \gamma', 'FontSize', fs);
ylabel('Stress', 'FontSize', fs);
title('Stress-Strain Curve and Yield Point', 'FontSize', fs + 2, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', fs - 2);
set(gca, 'FontSize', fs, 'LineWidth', 1.2);

% --- Panel 2: elastic energy curve (exact bond-sum vs Exp.7 numerical integration) ---
subplot(2,1,2); hold on; box on; grid on;
plot(gammaVec, elasticEnergy, '-^', 'LineWidth', lw, 'Color', [0.47 0.67 0.19], ...
    'MarkerFaceColor', [0.47 0.67 0.19], 'MarkerSize', 5, ...
    'DisplayName', 'exact, E = \Sigma 1/2 k(d_{new}-d_0)^2');
plot(gammaVec, elasticEnergy_numeric, '--x', 'LineWidth', lw, 'Color', [0.49 0.18 0.56], ...
    'MarkerSize', 7, 'DisplayName', 'numeric, \int\tau_{xy}\,d\gamma\cdotV (Exp.7)');
if ~isempty(gammaYieldRefined)
    xline(gammaYieldRefined, '-', 'LineWidth', lw, 'Color', [0.6 0 0], 'HandleVisibility', 'off');
end
xlabel('Engineering shear strain \gamma', 'FontSize', fs);
ylabel('Elastic energy', 'FontSize', fs);
title('Elastic Strain Energy: Exact Bond-Sum vs. Numerical Integration', 'FontSize', fs + 2, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', fs - 2);
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
