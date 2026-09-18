function modeShape = extractModeShape(eigenvectorColumn, mass, Natoms)
%EXTRACTMODESHAPE Turn a mass-weighted eigenvector of Dmat into real atomic displacements u = M^-1/2 v.
%  modeShape = extractModeShape(eigenvectorColumn, mass, Natoms)   mass: scalar or [Natoms x 1].
%  Returns [Natoms x 3] (x,y,z per atom), scaled so the largest atom displacement has length 1.
%  Eigenvector layout is [atom1_x; atom1_y; atom1_z; atom2_x; ...] as assembled in buildDynamicalMatrix.

    v = eigenvectorColumn(:);

    if length(v) ~= 3 * Natoms
        error('extractModeShape:badSize', ...
            'Eigenvector length (%d) does not match 3 * Natoms (%d).', length(v), 3 * Natoms);
    end

    modeShape = reshape(v, 3, Natoms)';

    if nargin >= 2 && ~isempty(mass)
        m = mass(:);
        if numel(m) == 1
            modeShape = modeShape / sqrt(m);
        elseif numel(m) == Natoms
            modeShape = bsxfun(@rdivide, modeShape, sqrt(m));
        end
    end

    magnitudes = sqrt(sum(modeShape.^2, 2));
    maxDisp = max(magnitudes);
    if maxDisp > 1e-12
        modeShape = modeShape / maxDisp;
    end
end
