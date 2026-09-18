function [latticeVectors, basisFrac, basisMasses] = parseCrystalFile(filepath)
%PARSECRYSTALFILE Read a simple text crystal definition (.txt/.dat) as an alternative to a .cif.
%  [latticeVectors, basisFrac, basisMasses] = parseCrystalFile(filepath)
%  Blank lines and lines starting with % or # are ignored.
%  First three numeric rows = Cartesian lattice vectors; every later row = "fracX fracY fracZ mass".

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
