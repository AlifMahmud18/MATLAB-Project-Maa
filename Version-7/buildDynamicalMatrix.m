function [K, M, Dmat] = buildDynamicalMatrix(pos, bonds, springConstant, masses)
%BUILDDYNAMICALMATRIX Assemble global Stiffness (K), Mass (M), and
%   Dynamical (Dmat) matrices for the mass-spring crystal network.
%
%   [K, M, Dmat] = buildDynamicalMatrix(pos, bonds, springConstant, masses)
%
%   PHYSICS BACKGROUND
%   Each atom has 3 degrees of freedom (x,y,z displacement from its
%   equilibrium lattice site). We model each bond as a CENTRAL-FORCE
%   spring: the restoring force only acts along the bond direction (this
%   is the harmonic approximation of a pairwise potential like a
%   Lennard-Jones bond, expanded to 2nd order).
%
%   For a bond between atoms i and j with unit vector u (pointing i->j)
%   and spring constant k, the potential energy for small displacements
%   (dr_i, dr_j) is:
%
%       U_bond = (1/2) * k * [u . (dr_i - dr_j)]^2
%
%   Expanding this gives 3x3 coupling blocks of the form k*(u*u'):
%
%       K(i,i) += k*(u*u')      K(j,j) += k*(u*u')
%       K(i,j) -= k*(u*u')      K(j,i) -= k*(u*u')
%
%   Stacking these blocks for every bond gives the global 3N x 3N
%   stiffness matrix K. The mass matrix M is diag(m1,m1,m1,m2,m2,m2,...)
%   - each atom's own mass repeated for its 3 DOF (x,y,z) - so different
%   atoms in the basis (e.g. a multi-species crystal) can carry different
%   masses.
%
%   The DYNAMICAL MATRIX is the mass-weighted stiffness matrix:
%
%       Dmat = M^(-1/2) * K * M^(-1/2)
%
%   Its eigenvalues are the squared angular frequencies (omega^2) of the
%   normal modes, and its eigenvectors, once "un-mass-weighted", give the
%   3D displacement pattern of each atom for that mode.
%
%   INPUTS
%     pos            - [Natoms x 3] atom positions (from generateLattice)
%     bonds          - [Nbonds x 5] bond list (from generateLattice)
%     springConstant - scalar spring constant k (same for every bond)
%     masses         - [Natoms x 1] atomic mass of each atom, aligned
%                       row-for-row with pos (e.g. from
%                       generateLatticeGeneral)
%
%   OUTPUTS
%     K    - [3N x 3N] global stiffness matrix
%     M    - [3N x 3N] diagonal mass matrix
%     Dmat - [3N x 3N] dynamical matrix (symmetric)

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
        u = bonds(b, 4:6);          % unit vector, already computed
        block = springConstant * (u' * u);   % 3x3 outer product

        idx_i = 3*(i-1) + (1:3);
        idx_j = 3*(j-1) + (1:3);

        K(idx_i, idx_i) = K(idx_i, idx_i) + block;
        K(idx_j, idx_j) = K(idx_j, idx_j) + block;
        K(idx_i, idx_j) = K(idx_i, idx_j) - block;
        K(idx_j, idx_i) = K(idx_j, idx_i) - block;
    end

    % Expand per-atom masses to 3 DOF each (x,y,z), in the same
    % atom-index order used above for idx_i/idx_j, so massVec(3*(i-1)+1:3)
    % all equal masses(i).
    massVec = kron(masses(:), [1; 1; 1]);
    M = diag(massVec);

    invSqrtM = diag(1 ./ sqrt(massVec));
    Dmat = invSqrtM * K * invSqrtM;

    % Numerical cleanup: force exact symmetry (roundoff can leave tiny
    % asymmetries that confuse the eigensolver later)
    Dmat = (Dmat + Dmat') / 2;
end