function [latticeVectors, basisFrac, basisMasses] = parseCrystalFile(filepath)
%PARSECRYSTALFILE Read a simple custom crystal-definition text file
%   (.txt/.dat), as an alternative to loading a full .cif via
%   parseCIFFile.m.
%
%   [latticeVectors, basisFrac, basisMasses] = parseCrystalFile(filepath)
%
%   FILE FORMAT
%   Blank lines and lines starting with '%' or '#' are ignored. The
%   first three numeric rows found are the lattice vectors (one 3D
%   Cartesian vector per row); every numeric row after that is one
%   basis atom as "fracX fracY fracZ mass". Example:
%
%       % Lattice vectors (Cartesian)
%       1.0 0.0 0.0
%       0.0 1.0 0.0
%       0.0 0.0 1.0
%       % Basis atoms: fracX fracY fracZ mass
%       0.0 0.0 0.0 1.0
%       0.5 0.5 0.5 2.0
%
%   INPUT
%     filepath - path to the .txt/.dat crystal file
%
%   OUTPUTS
%     latticeVectors - [3 x 3] Cartesian lattice vectors (rows)
%     basisFrac      - [Nbasis x 3] fractional coordinates
%     basisMasses    - [Nbasis x 1] atomic masses

    fid = fopen(filepath, 'r');
    if fid == -1
        error('parseCrystalFile:cannotOpen', 'Could not open file: %s', filepath);
    end

    latticeVectors = [];
    basisFrac = [];
    basisMasses = [];

    while true
        rawLine = fgetl(fid);
        if ~ischar(rawLine)
            break;
        end
        line = strtrim(rawLine);
        if isempty(line) || strncmp(line, '%', 1) || strncmp(line, '#', 1)
            continue;
        end

        values = str2double(strsplit(line));
        values = values(~isnan(values));

        if size(latticeVectors, 1) < 3 && numel(values) >= 3
            latticeVectors = [latticeVectors; values(1:3)]; %#ok<AGROW>
        elseif numel(values) >= 4
            basisFrac = [basisFrac; values(1:3)]; %#ok<AGROW>
            basisMasses = [basisMasses; values(4)]; %#ok<AGROW>
        end
    end
    fclose(fid);

    if size(latticeVectors, 1) ~= 3
        error('parseCrystalFile:badLattice', ...
            'Could not find 3 lattice-vector rows in %s.', filepath);
    end
    if isempty(basisFrac)
        error('parseCrystalFile:noAtoms', ...
            'Could not find any "fracX fracY fracZ mass" basis rows in %s.', filepath);
    end
end
