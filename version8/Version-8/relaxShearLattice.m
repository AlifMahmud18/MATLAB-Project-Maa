function [coords_relaxed, coords_affine, maxResidualBefore, maxResidualAfter, iterInfo] = ...
    relaxShearLattice(pos, bonds, springConstant, K_matrix, shearVal, solverMethod)
%RELAXSHEARLATTICE Equilibrium of the interior atoms when the top/bottom layers are sheared.
%  [coords_relaxed, coords_affine, resBefore, resAfter, iterInfo] = relaxShearLattice(pos, bonds, k, K, shearVal, 'gj')
%  Boundary layers (z_min, z_max) are fixed at the affine displacement u_x = gamma*(z - z_min).
%  Free atoms solve K_ff*u_f = -K_fb*u_b with gaussJordanSolve (EEE 212 Exp. 5).
%  resBefore/resAfter = max residual force on free atoms for the affine vs relaxed field.

    if nargin < 6 || isempty(solverMethod)
        solverMethod = 'gj';
    end

    Natoms = size(pos, 1);
    z = pos(:, 3);
    z_min = min(z);
    z_max = max(z);

    factor = (z - z_min);
    dispField_affine = zeros(Natoms, 3);
    dispField_affine(:, 1) = shearVal * factor;
    coords_affine = pos + dispField_affine;

    boundaryIdx = find(z == z_min | z == z_max);
    freeIdx = setdiff((1:Natoms)', boundaryIdx);

    dofOf = @(idx) reshape([3*(idx-1)+1, 3*(idx-1)+2, 3*(idx-1)+3]', [], 1);
    freeDofs = dofOf(freeIdx);
    boundaryDofs = dofOf(boundaryIdx);

    u_b = reshape(dispField_affine(boundaryIdx, :)', [], 1);

    gamma_eff = shearVal;
    S_xz = [0 0 1; 0 0 0; 0 0 0];
    [~, ~, ~, ~, ~, atomForce_affine] = computeShearForces(pos, bonds, springConstant, gamma_eff, S_xz);
    maxResidualBefore = max(max(abs(atomForce_affine(freeIdx, :))));

    K_ff = K_matrix(freeDofs, freeDofs);
    K_fb = K_matrix(freeDofs, boundaryDofs);
    rhs = -K_fb * u_b;

    iterInfo = struct('method', solverMethod);
    switch lower(solverMethod)
        case 'gj'
            u_f = gaussJordanSolve(K_ff, rhs);
        otherwise
            error('relaxShearLattice:badMethod', 'solverMethod must be ''gj''.');
    end

    dispField_relaxed = dispField_affine;
    dispField_relaxed(freeIdx, :) = reshape(u_f, 3, numel(freeIdx))';
    coords_relaxed = pos + dispField_relaxed;

    u_relaxed_full = reshape(dispField_relaxed', [], 1);
    forceAfter = K_matrix(freeDofs, :) * u_relaxed_full;
    maxResidualAfter = max(abs(forceAfter));
end
