%% runPhononShearAnalysis.m  (Version 7)
% Command-line version of the two app buttons: builds a crystal (built-in or
% ANY .cif), runs the shear/phonon sweeps, draws the 6 phonon-softening
% figures (P1-P6) and 2 bond-breaking figures (Q1-Q2), saves PNGs to
% shear_figures/, and prints validation checks.
% Set cifFile = '' to use a built-in crystal instead.

clear; clc; close all;
here = fileparts(mfilename('fullpath'));  addpath(here);
outDir = fullfile(here, 'shear_figures');

cifFile   = fullfile(here, 'CIF files', 'AL2O3.cif');
builtin   = 'bcc';
N         = 1;
% k = 1 is a reduced-unit spring constant: use a force constant fitted to the real material (e.g. from its elastic constants).
k         = 1;
gammaMax  = 0.5;
runAllCIFs = true;

[pos, bonds, masses, LV, basis] = buildLattice(cifFile, builtin, N, k);
fprintf('%d atoms, %d bonds\n', size(pos, 1), size(bonds, 1));

S = phononShearSuite(pos, bonds, masses, LV, basis, k, gammaMax);
plotShearSoftening(S, outDir);
plotShearBonds(S, outDir);
M = S.cases{S.main};
fprintf('first bond breaks at gamma = %.4f, lowest mode imaginary at gamma = %.4f\n', min(M.gammaFirstBreak), M.gammaInstab);
fprintf('max non-affine relaxation residual (force): %.2e\n', max(M.relaxResidual));

A = S.cases{1};
[~, sc] = min(abs(A.gamma - 0.30));  gChk = A.gamma(sc);
[K0, ~, ~] = buildDynamicalMatrix(pos, bonds, k, masses);
[~, ~, ~, freqsApp] = applyShearForce(pos, bonds, masses, K0, k, gChk);
wApp = sort(freqsApp(:));  wApp = wApp(4:end);
fprintf('Check 1: affine model vs applyShearForce.m at gamma=%.3f -> max |d omega| = %.2e\n', gChk, max(abs(sqrt(max(A.lambda(:, sc), 0)) - wApp)));
[~, ~, bv, ~, bfv] = computeShearForces(pos, bonds, k, gChk, [0 0 1; 0 0 0; 0 0 0]);
Vol = size(pos, 1) * abs(det(LV)) / size(basis, 1);
fprintf('Check 2: tau_xz affine model = %.6g vs computeShearForces.m virial = %.6g\n', A.tau(sc), sum(bfv(:, 1) .* bv(:, 3)) / Vol);
[~, ~, D0] = buildDynamicalMatrix(pos, bonds, k, masses);
w0app = sort(sqrt(max(eig(D0), 0)));
fprintf('Check 3: gamma=0 spectrum vs buildDynamicalMatrix (3 translations dropped) -> max |d omega| = %.2e\n', max(abs(w0app(4:end) - sqrt(max(A.lambda(:, 1), 0)))));

if runAllCIFs
    files = dir(fullfile(here, 'CIF files', '*.cif'));
    fprintf('\n%-12s %6s %6s %7s %10s %10s %12s\n', 'crystal', 'atoms', 'bonds', 'shells', 'gamma_break', 'gamma_inst', 'soft(last)/soft0');
    for f = files'
        try
            [p, b, m, lv, bs, sh] = buildLattice(fullfile(f.folder, f.name), '', 1, k);
            T = phononShearSuite(p, b, m, lv, bs, k, gammaMax, 30);
            R = T.cases{3};  nV = numel(R.gamma);
            if ~isnan(R.gammaInstab), nV = find(R.gamma < R.gammaInstab, 1, 'last'); end
            fprintf('%-12s %6d %6d %7d %10.3f %10.3f %12.3f\n', f.name, size(p, 1), size(b, 1), sh, min(R.gammaFirstBreak), R.gammaInstab, R.omegaSoft(nV) / R.omegaSoft(1));
        catch ME
            fprintf('%-12s failed: %s\n', f.name, ME.message);
        end
    end
end

function [pos, bonds, masses, LV, basis, shells] = buildLattice(cifFile, builtin, N, k)
    if ~isempty(cifFile)
        [LV, basis, bm] = parseCIFFile(cifFile);
    else
        [LV, basis, bm] = getBuiltinCrystalDef(builtin);
    end
    [pos, bonds, masses, shells] = autoRigidShells(LV, basis, bm, N, k);
end
