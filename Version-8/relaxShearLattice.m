function [coords_relaxed, coords_affine, maxResidualBefore, maxResidualAfter, iterInfo] = ...
    relaxShearLattice(pos, bonds, springConstant, K_matrix, shearVal, solverMethod)
%RELAXSHEARLATTICE Improve on the purely affine (Cauchy-Born) shear model
%   by solving for true mechanical equilibrium of the interior atoms,
%   using the linear-system methods of EEE 212 Exp. 5.
%
%   [coords_relaxed, coords_affine, maxResidualBefore, maxResidualAfter, iterInfo] = ...
%       relaxShearLattice(pos, bonds, springConstant, K_matrix, shearVal, solverMethod)
%
%   MOTIVATION (see the note in computeShearForces.m / applyShearForce.m):
%   applyShearForce.m imposes the SAME affine displacement
%       u_x(y) = shearVal * (y - y_min) / (y_max - y_min),  u_y = u_z = 0
%   on every atom, including interior ones. That is exact only for an
%   infinite, homogeneous crystal; in a finite supercell the free
%   surfaces in y leave interior atoms with a non-zero net Hooke's-law
%   force under the harmonic stiffness matrix K ("surface relaxation").
%
%   This function instead treats the boundary layers (y = y_min and
%   y = y_max, i.e. the atoms a shear rig would physically grip) as
%   DIRICHLET boundary conditions fixed at the affine displacement, and
%   solves for the equilibrium displacement of every other ("free")
%   atom by requiring zero net force on each of them:
%
%       K_ff * u_f = -K_fb * u_b        (partition boundary/free DOFs)
%
%   which is a straightforward n-unknown linear system - solved here
%   with gaussJordanSolve.m or gaussSeidelSolve.m (Exp. 5) instead of a
%   second affine assumption. This is well-posed regardless of periodic
%   wraparound: K_ff/K_fb are plain submatrices of the already-correct
%   global K, and u_b assigns each BOUNDARY atom its own prescribed
%   value directly (no relative/difference computation involved), so the
%   solved u_f is exactly the free-atom displacement that zeroes every
%   free atom's net bond force GIVEN those prescribed boundary values.
%
%   INPUTS
%     pos            - [Natoms x 3] pristine (unsheared) positions
%     bonds          - [Nbonds x 6] bond list (from generateLatticeGeneral.m)
%     springConstant - scalar spring constant k
%     K_matrix       - [3N x 3N] stiffness matrix at the pristine
%                       configuration (from buildDynamicalMatrix.m).
%                       Using the pristine K is the same small-
%                       displacement/harmonic linearization already used
%                       everywhere else in this codebase.
%     shearVal       - scalar shear magnitude (same convention as
%                       applyShearForce.m)
%     solverMethod   - (optional) 'gj' for gaussJordanSolve (default) or
%                       'gs' for gaussSeidelSolve
%
%   OUTPUTS
%     coords_relaxed     - [Natoms x 3] boundary atoms at the affine
%                           displacement, free atoms at the solved
%                           equilibrium displacement
%     coords_affine       - [Natoms x 3] the old fully-affine result,
%                           returned for side-by-side comparison
%     maxResidualBefore   - max |net force| over free atoms under the
%                           fully affine field (the problem being fixed)
%     maxResidualAfter    - max |net force| over free atoms after
%                           relaxation (should be ~0 by construction)
%     iterInfo            - struct with .method and, for Gauss-Seidel,
%                           .iters / .converged

    if nargin < 6 || isempty(solverMethod)
        solverMethod = 'gj';
    end

    Natoms = size(pos, 1);
    y = pos(:, 2);
    y_min = min(y);
    y_max = max(y);
    span_y = max(y_max - y_min, eps);

    factor = (y - y_min) / span_y;               % [Natoms x 1], 0 at bottom, 1 at top
    dispField_affine = zeros(Natoms, 3);
    dispField_affine(:, 1) = shearVal * factor;   % same profile as applyShearForce.m
    coords_affine = pos + dispField_affine;

    boundaryIdx = find(y == y_min | y == y_max);
    freeIdx = setdiff((1:Natoms)', boundaryIdx);

    dofOf = @(idx) reshape([3*(idx-1)+1, 3*(idx-1)+2, 3*(idx-1)+3]', [], 1);
    freeDofs = dofOf(freeIdx);
    boundaryDofs = dofOf(boundaryIdx);

    u_b = reshape(dispField_affine(boundaryIdx, :)', [], 1);

    % Residual net force BEFORE relaxing: apply the affine field to every
    % atom (including "free" ones) and see how far each free atom's net
    % bond force is from zero.
    %
    % This is deliberately computed via computeShearForces.m's own
    % (already wraparound-safe) atomForce output rather than
    % "K_matrix(freeDofs,:) * u_affine_full": generateLatticeGeneral.m's
    % bonds can connect an atom to a periodic IMAGE of its neighbor, and
    % a single global displacement array built from the raw affine
    % formula u=gamma*S*r does not correctly represent the relative
    % displacement of such a wrapped bond (same issue documented in
    % computeShearForces.m and fixed there for exactly this reason) -
    % reusing that already-correct machinery here sidesteps it.
    gamma_eff = shearVal / span_y;    % this file's dx=shearVal*(y-ymin)/span_y is gamma_eff*y up to an additive constant that cancels in every bond difference
    [~, ~, ~, ~, ~, atomForce_affine] = computeShearForces(pos, bonds, springConstant, gamma_eff);
    maxResidualBefore = max(max(abs(atomForce_affine(freeIdx, :))));

    K_ff = K_matrix(freeDofs, freeDofs);
    K_fb = K_matrix(freeDofs, boundaryDofs);
    rhs = -K_fb * u_b;

    iterInfo = struct('method', solverMethod);
    switch lower(solverMethod)
        case 'gj'
            u_f = gaussJordanSolve(K_ff, rhs);
        case 'gs'
            [u_f, iters, converged] = gaussSeidelSolve(K_ff, rhs);
            iterInfo.iters = iters;
            iterInfo.converged = converged;
        otherwise
            error('relaxShearLattice:badMethod', 'solverMethod must be ''gj'' or ''gs''.');
    end

    dispField_relaxed = dispField_affine;
    dispField_relaxed(freeIdx, :) = reshape(u_f, 3, numel(freeIdx))';
    coords_relaxed = pos + dispField_relaxed;

    % Residual net force AFTER relaxing: this is exactly the quantity
    % the linear solve above was built to zero out (K_ff*u_f + K_fb*u_b),
    % using the SAME index-based K_matrix*u operation the solve itself
    % used - so, unlike "before", this is self-consistent by
    % construction regardless of wraparound and should land at ~0.
    u_relaxed_full = reshape(dispField_relaxed', [], 1);
    forceAfter = K_matrix(freeDofs, :) * u_relaxed_full;
    maxResidualAfter = max(abs(forceAfter));
end
