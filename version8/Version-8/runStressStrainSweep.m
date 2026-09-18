% runStressStrainSweep.m - shear-strain sweep of a mass-spring lattice: stress, energy, shear modulus G, yield strain.
% Each strain: computeShearForces -> computeAtomicStress -> checkPlasticity (von Mises / Tresca).
% EEE 212 tools: least-squares fit + gaussJordanSolve (G), compositeTrapz (energy), bisection/Newton (yield strain),
% Newton divided-difference interpolation check; the four small routines are local functions at the bottom.
% Figure: stress-strain curve with yield marker, and exact vs numerically integrated elastic energy.

clear; clc; close all;

% ASSUMED reduced units: a0, k and mass are generic; use the real lattice constant, force constant and mass for a real material.
a0 = 1.0;
latticeVectors = a0 * eye(3);
basisFrac      = [0 0 0];
basisMasses    = 1.0;
N         = 4;
numShells = 2;

[pos, bonds, masses] = generateLatticeGeneral(latticeVectors, basisFrac, basisMasses, N, numShells);
Natoms = size(pos, 1);

springConstant = 1.0;

cellVolume   = abs(det(latticeVectors));
atomsPerCell = size(basisFrac, 1);
atomicVolume = cellVolume / atomsPerCell;
totalVolume  = cellVolume * N^3;

% ILLUSTRATIVE yield strengths: calibrate to the real material (e.g. literature values, Callister Ch. 6).
sigmaYield = 0.15;
tauYield   = 0.10;

gammaMax = 0.30;
nSteps   = 60;
gammaVec = linspace(0, gammaMax, nSteps)';
h_gamma  = gammaVec(2) - gammaVec(1);

avgVonMises   = zeros(nSteps, 1);
avgSigmaXY    = zeros(nSteps, 1);
elasticEnergy = zeros(nSteps, 1);
anyYieldVM    = false(nSteps, 1);
anyYieldTr    = false(nSteps, 1);

for s = 1:nSteps
    gamma = gammaVec(s);

    [~, ~, bondVecs_sh, bondForceMag, bondForceVec, ~] = ...
        computeShearForces(pos, bonds, springConstant, gamma);

    [stressTensor, principalStresses, vonMises] = ...
        computeAtomicStress(bonds, bondVecs_sh, bondForceVec, Natoms, atomicVolume);

    [vmYield, trYield] = checkPlasticity(principalStresses, vonMises, sigmaYield, tauYield);
    anyYieldVM(s) = any(vmYield);
    anyYieldTr(s) = any(trYield);

    avgVonMises(s) = mean(vonMises);
    avgSigmaXY(s)  = mean(squeeze(stressTensor(1, 2, :)));

    elasticEnergy(s) = 0.5 * sum(bondForceMag.^2) / springConstant;
end

firstYieldIdx = find(anyYieldVM | anyYieldTr, 1);
if isempty(firstYieldIdx)
    linearIdx = 1:nSteps;
else
    linearIdx = 1:max(firstYieldIdx - 1, 2);
end

cFit  = polyLeastSquaresFit(gammaVec(linearIdx), avgSigmaXY(linearIdx), 1);
G_fit = cFit(2);
fprintf('Fitted shear modulus G = %.6g (from tau_xy = G*gamma, %d-point linear fit, Exp.4/Exp.5)\n', ...
    G_fit, numel(linearIdx));

tau_ideal = G_fit / (2 * pi);
fprintf('Frenkel ideal shear strength estimate: tau_ideal ~= G/(2*pi) = %.6g\n', tau_ideal);

[~, energyDensityCum] = compositeTrapz(avgSigmaXY, h_gamma);
elasticEnergy_numeric = energyDensityCum * totalVolume;

energyCheckErr = max(abs(elasticEnergy_numeric - elasticEnergy));
energyCheckRel = energyCheckErr / max(abs(elasticEnergy(end)), eps);
fprintf(['Exp.7 cross-check: trapezoidal-integrated energy vs exact bond-sum energy -> ' ...
    'max abs diff = %.3e (relative %.3e)\n'], energyCheckErr, energyCheckRel);

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

lw = 2; fs = 13;

figure('Color', 'w', 'Position', [100 100 950 750]);

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

function [xMid, yMid, iters] = bisectionRoot(f, xLower, xUpper, xTol)
    if nargin < 4 || isempty(xTol), xTol = 1e-8; end
    yLower = f(xLower);
    yUpper = f(xUpper);
    if yLower * yUpper > 0
        error('bisectionRoot:noBracket', 'f(xLower) and f(xUpper) must have opposite signs.');
    end
    xMid = (xLower + xUpper) / 2.0;
    yMid = f(xMid);
    iters = 0;
    while (xUpper - xLower) / 2.0 > xTol
        iters = iters + 1;
        if yLower * yMid > 0
            xLower = xMid;
            yLower = yMid;
        else
            xUpper = xMid;
        end
        xMid = (xLower + xUpper) / 2.0;
        yMid = f(xMid);
    end
end

function [root, iters] = newtonRaphsonRoot(f, x0, xTol, fDeriv)
    if nargin < 3 || isempty(xTol), xTol = 1e-8; end
    if nargin < 4 || isempty(fDeriv), fDeriv = @(x) centralDiffDerivative(f, x); end
    iters = 1;
    dx = -f(x0) / fDeriv(x0);
    root = x0 + dx;
    maxIter = 200;
    while abs(dx) > xTol && iters < maxIter
        dx = -f(root) / fDeriv(root);
        root = root + dx;
        iters = iters + 1;
    end
end

function yq = newtonDividedDiffInterp(x, y, xq)
    x = x(:);
    y = y(:);
    n = numel(x);
    if numel(y) ~= n
        error('newtonDividedDiffInterp:sizeMismatch', 'x and y must have the same length.');
    end
    D = zeros(n, n);
    D(:, 1) = y;
    for j = 2:n
        for k = j:n
            D(k, j) = (D(k, j-1) - D(k-1, j-1)) / (x(k) - x(k-j+1));
        end
    end
    coeffs = diag(D);
    yq = coeffs(n) * ones(size(xq));
    for k = n-1:-1:1
        yq = coeffs(k) + (xq - x(k)) .* yq;
    end
end

function coeffs = polyLeastSquaresFit(x, y, degree)
    x = x(:);
    y = y(:);
    d = degree;
    powerSums = zeros(1, 2*d + 1);
    for p = 0:2*d
        powerSums(p+1) = sum(x.^p);
    end
    M = zeros(d+1, d+1);
    rhs = zeros(d+1, 1);
    for j = 0:d
        for k = 0:d
            M(j+1, k+1) = powerSums(j+k+1);
        end
        rhs(j+1) = sum((x.^j) .* y);
    end
    coeffs = gaussJordanSolve(M, rhs);
end
