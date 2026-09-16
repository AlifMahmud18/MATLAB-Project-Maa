function [pos, bonds, masses] = generateLatticeGeneral(varargin)
%GENERATELATTICEGENERAL Build atom positions, bonds with PBC, and masses.
%
% Usage:
%   [pos, bonds, masses] = generateLatticeGeneral(cifFilename, N, numShells)
%   [pos, bonds, masses] = generateLatticeGeneral(latticeVectors, basisFrac, basisMasses, N, numShells)

    % --- Input Parsing & Polymorphism ---
    if ischar(varargin{1}) || isstring(varargin{1})
        % Called via CIF File Path: (cifFilename, [N], [numShells])
        cifFilename = varargin{1};
        N = 1;
        numShells = 1;
        if nargin >= 2 && ~isempty(varargin{2})
            N = varargin{2};
        end
        if nargin >= 3 && ~isempty(varargin{3})
            numShells = varargin{3};
        end
        [latticeVectors, basisFrac, basisMasses] = parseCIFFile(cifFilename);
    else
        % Called via Direct Matrices: (latticeVectors, basisFrac, basisMasses, [N], [numShells])
        latticeVectors = varargin{1};
        basisFrac = varargin{2};
        basisMasses = varargin{3};
        N = 1;
        numShells = 1;
        if nargin >= 4 && ~isempty(varargin{4})
            N = varargin{4};
        end
        if nargin >= 5 && ~isempty(varargin{5})
            numShells = varargin{5};
        end
    end

    Mbasis = size(basisFrac, 1);

    % --- Tile basis across N x N x N supercell ---
    superLattice = N * latticeVectors; % 3x3 supercell lattice matrix
    pos = [];
    for ix = 0:N-1
        for iy = 0:N-1
            for iz = 0:N-1
                cellOriginCart = [ix, iy, iz] * latticeVectors;
                cellOriginExpanded = repmat(cellOriginCart, Mbasis, 1);
                basisCart = basisFrac * latticeVectors;
                pos = [pos; cellOriginExpanded + basisCart]; %#ok<AGROW>
            end
        end
    end

    Natoms = size(pos, 1);
    nCells = N^3;
    masses = repmat(basisMasses(:), nCells, 1);

    % --- Pairwise distances using Minimum Image Convention (PBC) ---
    fracPos = pos / superLattice; % Supercell fractional coordinates
    D = zeros(Natoms, Natoms);
    
    for i = 1:Natoms
        df = fracPos - repmat(fracPos(i, :), Natoms, 1);
        df = df - round(df); % Wrap fractional displacement to [-0.5, 0.5]
        dr = df * superLattice; % Minimum image displacement vector in Cartesian space
        D(:, i) = sqrt(sum(dr.^2, 2));
    end
    
    D(D < 1e-8) = Inf; % Ignore self-distances

    % --- Identify neighbor distance shells ---
    roundedD = round(D(:) * 1e4) / 1e4;
    shellDistances = unique(roundedD(isfinite(roundedD)));
    shellDistances = sort(shellDistances);
    if isempty(shellDistances)
        error('generateLatticeGeneral:noNeighbors', 'No neighbor interactions found in lattice.');
    end
    if numShells > length(shellDistances)
        numShells = length(shellDistances);
    end
    cutoff = shellDistances(numShells) * 1.001;

    % --- Build bond network under Periodic Boundary Conditions ---
    bonds = [];
    for i = 1:Natoms
        for j = i+1:Natoms
            df = fracPos(j, :) - fracPos(i, :);
            df = df - round(df); % MIC wrapping
            dr = df * superLattice;
            dij = norm(dr);
            
            if dij <= cutoff
                uvec = dr / dij;
                bonds = [bonds; i, j, dij, uvec]; %#ok<AGROW>
            end
        end
    end
end