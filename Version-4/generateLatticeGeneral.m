function [pos, bonds, masses] = generateLatticeGeneral(latticeVectors, basisFrac, basisMasses, N, numShells)
%GENERATELATTICEGENERAL Build atom positions, bonds, and masses for a crystal supercell.

    % --- Defaults & Validation ---
    if nargin < 5 || isempty(numShells)
        numShells = 1;
    end

    if ~isequal(size(latticeVectors), [3 3])
        error('generateLatticeGeneral:badLatticeVectors', ...
            'latticeVectors must be a 3x3 matrix (rows = a1, a2, a3).');
    end

    if numel(basisMasses) ~= size(basisFrac, 1)
        error('generateLatticeGeneral:badBasisMasses', ...
            'basisMasses must have one entry per row of basisFrac (%d expected, got %d).', ...
            size(basisFrac, 1), numel(basisMasses));
    end

    % Convert basis coordinates to Cartesian
    basisCart = basisFrac * latticeVectors;   % Mx3
    Mbasis = size(basisCart, 1);

    % 1. Tile basis across N x N x N grid of unit cells
    pos = [];
    for ix = 0:N-1
        for iy = 0:N-1
            for iz = 0:N-1
                cellOriginCart = [ix, iy, iz] * latticeVectors;
                cellOriginExpanded = repmat(cellOriginCart, Mbasis, 1);
                pos = [pos; cellOriginExpanded + basisCart]; %#ok<AGROW>
            end
        end
    end

    Natoms = size(pos, 1);

    % 2. Replicate basis masses so masses array matches Natoms row-for-row
    nCells = N^3;
    masses = repmat(basisMasses(:), nCells, 1);

    % 3. Pairwise distances
    D = zeros(Natoms, Natoms);
    for i = 1:Natoms
        diffs = pos - repmat(pos(i, :), Natoms, 1);
        D(:, i) = sqrt(sum(diffs.^2, 2));
    end
    D(D < 1e-8) = Inf;   % ignore self-distances

    % 4. Identify neighbor distance shells
    roundedD = round(D(:) * 1e4) / 1e4;
    shellDistances = unique(roundedD(isfinite(roundedD)));
    shellDistances = sort(shellDistances);

    if numShells > length(shellDistances)
        numShells = length(shellDistances);
    end
    cutoff = shellDistances(numShells) * 1.001;

    % 5. Build bond network
    bonds = [];
    for i = 1:Natoms
        for j = i+1:Natoms
            dij = D(i, j);
            if dij <= cutoff
                uvec = (pos(j, :) - pos(i, :)) / dij;
                bonds = [bonds; i, j, dij, uvec]; %#ok<AGROW>
            end
        end
    end
end