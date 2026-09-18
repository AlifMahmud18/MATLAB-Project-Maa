function [K, M, Dmat] = buildDynamicalMatrix(pos, bonds, springConstant, masses)
%BUILDDYNAMICALMATRIX Stiffness K, mass M and dynamical matrix Dmat = M^-1/2 K M^-1/2 of the crystal.
%  [K, M, Dmat] = buildDynamicalMatrix(pos, bonds, springConstant, masses)
%  Each bond is a central-force spring: block k*u*u' (u = bond unit vector) added to (i,i),(j,j), subtracted from (i,j),(j,i).
%  masses has one value per atom; the eigenvalues of Dmat are omega^2 of the normal modes.
%  Dmat is symmetrised at the end to remove round-off asymmetry.

    Natoms = size(pos, 1);

    if numel(masses) ~= Natoms
        error('buildDynamicalMatrix:badMasses', ...
            'masses must have one entry per atom (%d expected, got %d).', ...
            Natoms, numel(masses));
    end

    K = zeros(3 * Natoms);

    for b = 1:size(bonds, 1)
        i = bonds(b, 1);
        j = bonds(b, 2);
        u = bonds(b, 4:6);
        block = springConstant * (u' * u);

        idx_i = 3*(i-1) + (1:3);
        idx_j = 3*(j-1) + (1:3);

        K(idx_i, idx_i) = K(idx_i, idx_i) + block;
        K(idx_j, idx_j) = K(idx_j, idx_j) + block;
        K(idx_i, idx_j) = K(idx_i, idx_j) - block;
        K(idx_j, idx_i) = K(idx_j, idx_i) - block;
    end

    massVec = kron(masses(:), [1; 1; 1]);
    M = diag(massVec);

    invSqrtM = diag(1 ./ sqrt(massVec));
    Dmat = invSqrtM * K * invSqrtM;

    Dmat = (Dmat + Dmat') / 2;
end
