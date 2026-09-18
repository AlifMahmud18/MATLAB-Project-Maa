function [pos, bonds, masses] = generateLatticeGeneral(varargin)
%GENERATELATTICEGENERAL Tile the unit cell into an N x N x N supercell and build the periodic spring bonds.
%  [pos, bonds, masses] = generateLatticeGeneral(cifFilename, N, numShells)
%  [pos, bonds, masses] = generateLatticeGeneral(latticeVectors, basisFrac, basisMasses, N, numShells)
%  bonds rows are [i j d0 ux uy uz]: every periodic image within the first numShells neighbour distances.
%  Neighbours come from findPeriodicNeighborCandidates (Gamma-point model: all images move identically).

    if ischar(varargin{1}) || isstring(varargin{1})
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

    superLattice = N * latticeVectors;
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

    [candI, candJ, candDist, candU, shellDistances] = ...
        findPeriodicNeighborCandidates(pos, superLattice);

    if isempty(shellDistances)
        error('generateLatticeGeneral:noNeighbors', 'No neighbor interactions found in lattice.');
    end
    if numShells > length(shellDistances)
        numShells = length(shellDistances);
    end
    cutoff = shellDistances(numShells) * 1.001;

    withinCutoff = candDist <= cutoff;
    bonds = [candI(withinCutoff), candJ(withinCutoff), candDist(withinCutoff), candU(withinCutoff, :)];
end
