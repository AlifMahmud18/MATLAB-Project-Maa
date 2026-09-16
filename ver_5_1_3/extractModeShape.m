function modeShape = extractModeShape(eigenvectorColumn, mass, Natoms)
%EXTRACTMODESHAPE Convert one eigenvector of the DYNAMICAL matrix into a
%   physical atomic displacement pattern, ready for animation.
%
%   modeShape = extractModeShape(eigenvectorColumn, mass, Natoms)
%
%   WHY THIS CONVERSION IS NEEDED
%   Remember from buildDynamicalMatrix.m: Dmat = M^(-1/2) * K * M^(-1/2).
%   Its eigenvectors are in "mass-weighted" coordinates, not real physical
%   displacements. To get the actual displacement pattern u, you have to
%   undo the mass-weighting:
%
%       u = M^(-1/2) * v         (v = eigenvector of Dmat)
%
%   INPUTS
%     eigenvectorColumn - [3N x 1] one column of V from jacobiEigenSolver
%     mass              - scalar mass OR [Natoms x 1] vector of atomic masses
%     Natoms            - number of atoms (3N = length of eigenvectorColumn)
%
%   OUTPUT
%     modeShape - [Natoms x 3] displacement DIRECTION for each atom,
%                 normalized so the largest atomic displacement has
%                 magnitude 1.
%
%   NOTE ON ORDERING: eigenvectorColumn is laid out as
%   [atom1_x; atom1_y; atom1_z; atom2_x; atom2_y; atom2_z; ...] because
%   that's how buildDynamicalMatrix.m assembled the indices
%   (idx_i = 3*(i-1) + (1:3)). reshape(...,3,Natoms)' unpacks that
%   correctly back into one row per atom.

    v = eigenvectorColumn(:); % Force column vector

    if length(v) ~= 3 * Natoms
        error('extractModeShape:badSize', ...
            'Eigenvector length (%d) does not match 3 * Natoms (%d).', length(v), 3 * Natoms);
    end

    % 1. Unpack 3N x 1 vector into Natoms x 3 matrix (x, y, z per atom)
    modeShape = reshape(v, 3, Natoms)';   % [Natoms x 3]

    % 2. Undo mass-weighting (handles both scalar mass and per-atom mass vector)
    if nargin >= 2 && ~isempty(mass)
        m = mass(:);
        if numel(m) == 1
            modeShape = modeShape / sqrt(m);
        elseif numel(m) == Natoms
            % Use bsxfun instead of a bare "./" so this broadcasts an
            % [Natoms x 1] mass vector across the 3 (x,y,z) columns even
            % on MATLAB releases without implicit expansion (pre-R2016b) -
            % without it this throws "Matrix dimensions must agree".
            modeShape = bsxfun(@rdivide, modeShape, sqrt(m));
        end
    end

    % 3. Normalize so the largest single-atom displacement has length 1.
    magnitudes = sqrt(sum(modeShape.^2, 2));
    maxDisp = max(magnitudes);
    if maxDisp > 1e-12
        modeShape = modeShape / maxDisp;
    end
end