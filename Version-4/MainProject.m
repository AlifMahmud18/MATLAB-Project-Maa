%% MainProject.m
% EEE 212 Project: Interactive MATLAB-Based Crystalline Structure Viewer
% and Vibrational Mode Animator
%
% This is the ONE script to run. It calls out to helper FUNCTION files
% (generateLattice.m, buildDynamicalMatrix.m, jacobiEigenSolver.m,
% extractModeShape.m, animateVibration.m, plotLattice.m) which must all
% sit in the same folder as this script - MATLAB finds them by filename,
% the same way your Exp 1&2 lab sheet's main file called cal_pow.m.
%
% Run this whole file (F5). It will:
%   1. Build the crystal lattice and show a static 3D picture
%   2. Assemble the dynamical matrix
%   3. Solve it with the custom Jacobi eigensolver
%   4. Validate the solver against MATLAB's built-in eig()
%   5. Check the zero-mode count (Maxwell rigidity sanity check)
%   6. Animate a chosen vibrational mode

clear; clc; close all;

%% ---------------------------------------------------------------------
%  SECTION 1: Parameters - change these to explore the model
%  ------------------------------------------------------------------
structureType  = 'bcc';   % 'sc', 'bcc', or 'fcc'
N              = 2;       % NxNxN conventional cells
a              = 1.0;     % lattice constant (arbitrary units)
springConstant = 1.0;     % spring constant k (arbitrary units)
mass           = 1.0;     % atomic mass (arbitrary units)

% Bond network: how many neighbor-distance shells to connect with
% springs. 2 includes 1st + 2nd nearest neighbors, which keeps the
% network properly rigid (see Section 5 for why this matters). This
% replaced an earlier fixed-multiplier approach that turned out not to
% generalize correctly to every crystal structure (e.g. diamond) - see
% the comments in generateLatticeGeneral.m for the full explanation.
numShells = 2;

% Which vibrational mode to animate at the end (modes 1-6 are always
% the trivial zero-frequency rigid-body translations/rotations, so the
% first genuinely vibrating mode is usually index 7).
modeIndexToAnimate = 7;

% Animation appearance
amplitude = 0.25;
duration  = 8;      % seconds
fps       = 30;

%% ---------------------------------------------------------------------
%  SECTION 2: Build the lattice and show the static structure
%  ------------------------------------------------------------------
[pos, bonds, masses] = generateLattice(structureType, N, a, numShells);
Natoms = size(pos, 1);

figure('Name', 'Crystal Lattice (static)');
plotLattice(pos, bonds);
title(sprintf('%s lattice, %dx%dx%d cells, %d atoms', ...
    upper(structureType), N, N, N, Natoms));

%% ---------------------------------------------------------------------
%  SECTION 3: Build the dynamical matrix
%  ------------------------------------------------------------------
[K, M, Dmat] = buildDynamicalMatrix(pos, bonds, springConstant, masses);
fprintf('\nDynamical matrix size: %d x %d\n', size(Dmat,1), size(Dmat,2));

%% ---------------------------------------------------------------------
%  SECTION 4: Solve with the custom Jacobi eigensolver, validate vs eig()
%  ------------------------------------------------------------------
tic;
[V, eigVals] = jacobiEigenSolver(Dmat, 1e-12, 200);
tCustom = toc;

[eigVals, order] = sort(eigVals);
V = V(:, order);
fprintf('Custom Jacobi solver: %.4f sec\n', tCustom);

tic;
[~, D_ref] = eig(Dmat);
tRef = toc;
eigVals_ref = sort(diag(D_ref));

maxError = max(abs(eigVals - eigVals_ref));
fprintf('Built-in eig():       %.4f sec\n', tRef);
fprintf('\n>>> Max eigenvalue error (custom vs. built-in): %.3e <<<\n', maxError);

if maxError < 1e-6
    fprintf('VALIDATION PASSED: your Jacobi solver agrees with eig().\n');
else
    fprintf('VALIDATION FAILED: check your jacobiEigenSolver implementation.\n');
end

%% ---------------------------------------------------------------------
%  SECTION 5: Physical sanity check - zero modes (Maxwell rigidity)
%  ------------------------------------------------------------------
% A free-standing central-force cluster always has >= 6 zero-frequency
% modes (3 translations + 3 rotations - moving/rotating the WHOLE
% cluster rigidly stretches no bond). With numShells = 2 we included
% 1st+2nd neighbor bonds specifically so the network is fully rigid and
% this count comes out to EXACTLY 6, not more. (Try changing numShells
% to 1 - nearest neighbors only - and re-running: you'll see extra
% "floppy modes" from under-coordinated surface atoms.)

numZeroModes = sum(eigVals < 1e-6);
avgCoordination = 2 * size(bonds,1) / Natoms;
fprintf('\nAverage coordination number: %.2f\n', avgCoordination);
fprintf('Zero-frequency modes found: %d (should be exactly 6)\n', numZeroModes);

omega = sqrt(max(eigVals, 0));   % eigenvalues of Dmat = omega^2

figure('Name', 'Vibrational Frequency Spectrum');
stem(omega, 'filled');
xlabel('Mode index'); ylabel('\omega (angular frequency)');
title('Vibrational normal mode frequencies');
grid on;

fprintf('\nFirst 8 frequencies (omega):\n');
disp(omega(1:8)');

%% ---------------------------------------------------------------------
%  SECTION 6: Animate the chosen mode
%  ------------------------------------------------------------------
fprintf('\nAnimating mode %d, omega = %.4f\n', ...
    modeIndexToAnimate, omega(modeIndexToAnimate));

modeShape = extractModeShape(V(:, modeIndexToAnimate), masses, Natoms);

animateVibration(pos, bonds, modeShape, omega(modeIndexToAnimate), ...
    amplitude, duration, fps);

fprintf('\nDone. To see a different mode, change modeIndexToAnimate in\n');
fprintf('Section 1 above and re-run this whole script.\n');
