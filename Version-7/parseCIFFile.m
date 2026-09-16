function [latticeVectors, basisFrac, basisMasses] = parseCIFFile(filepath)
%PARSECIFFILE Read a Crystallographic Information File (.cif) with symmetry support.
%
%   [latticeVectors, basisFrac, basisMasses] = parseCIFFile(filepath)

    fid = fopen(filepath, 'r');
    if fid == -1
        error('parseCIFFile:cannotOpen', 'Could not open file: %s', filepath);
    end

    a = []; b = []; c = [];
    alpha = []; beta = []; gamma = [];
    asymmetricBasis = [];
    asymmetricMasses = [];
    symOps = {};

    inAtomLoop = false;
    inSymLoop = false;
    xCol = 0; yCol = 0; zCol = 0; symCol = 0; elemCol = 0;
    nAtomLoopCols = 0; nSymLoopCols = 0;

    while true
        rawLine = fgetl(fid);
        if ~ischar(rawLine)
            break;
        end
        line = strtrim(rawLine);
        if isempty(line) || strncmp(line, '#', 1)
            inAtomLoop = false;
            inSymLoop = false;
            continue;
        end

        lowerLine = lower(line);

        % --- Cell parameters ---
        if strncmpi(lowerLine, '_cell_length_a', 14)
            a = readCifNumber(line);
        elseif strncmpi(lowerLine, '_cell_length_b', 14)
            b = readCifNumber(line);
        elseif strncmpi(lowerLine, '_cell_length_c', 14)
            c = readCifNumber(line);
        elseif strncmpi(lowerLine, '_cell_angle_alpha', 17)
            alpha = readCifNumber(line);
        elseif strncmpi(lowerLine, '_cell_angle_beta', 16)
            beta = readCifNumber(line);
        elseif strncmpi(lowerLine, '_cell_angle_gamma', 17)
            gamma = readCifNumber(line);

        % --- Start of a loop_ block ---
        elseif strcmpi(line, 'loop_')
            inAtomLoop = false;
            inSymLoop = false;
            loopTags = {};
            headerPos = ftell(fid);

            while true
                tagLine = fgetl(fid);
                if ~ischar(tagLine)
                    break;
                end
                tagLineTrim = strtrim(tagLine);
                if isempty(tagLineTrim)
                    continue;
                end
                if strncmp(tagLineTrim, '_', 1)
                    loopTags{end+1} = lower(tagLineTrim); %#ok<AGROW>
                    headerPos = ftell(fid);
                else
                    fseek(fid, headerPos, 'bof');   % rewind to first data row
                    break;
                end
            end

            % Check for Atom site coordinates
            xIdx = find(strcmp(loopTags, '_atom_site_fract_x'), 1);
            yIdx = find(strcmp(loopTags, '_atom_site_fract_y'), 1);
            zIdx = find(strcmp(loopTags, '_atom_site_fract_z'), 1);

            % Check for Symmetry operations
            symIdx = find(strcmp(loopTags, '_space_group_symop_operation_xyz') | ...
                          strcmp(loopTags, '_symmetry_equiv_pos_as_xyz'), 1);

            % Check for atom type/element column (prefer explicit type symbol)
            typeIdx  = find(strcmp(loopTags, '_atom_site_type_symbol'), 1);
            labelIdx = find(strcmp(loopTags, '_atom_site_label'), 1);

            if ~isempty(xIdx) && ~isempty(yIdx) && ~isempty(zIdx)
                inAtomLoop = true;
                xCol = xIdx; yCol = yIdx; zCol = zIdx;
                if ~isempty(typeIdx)
                    elemCol = typeIdx;
                elseif ~isempty(labelIdx)
                    elemCol = labelIdx;
                else
                    elemCol = 0;
                end
                nAtomLoopCols = numel(loopTags);
            elseif ~isempty(symIdx)
                inSymLoop = true;
                symCol = symIdx;
                nSymLoopCols = numel(loopTags);
            end

        % --- Data row inside confirmed atom_site loop ---
        elseif inAtomLoop
            tokens = strsplit(line);
            tokens = tokens(~cellfun(@isempty, tokens));
            if numel(tokens) < nAtomLoopCols
                inAtomLoop = false;
                continue;
            end
            fx = str2double(stripUncertainty(tokens{xCol}));
            fy = str2double(stripUncertainty(tokens{yCol}));
            fz = str2double(stripUncertainty(tokens{zCol}));
            if any(isnan([fx fy fz]))
                inAtomLoop = false;
                continue;
            end
            asymmetricBasis = [asymmetricBasis; fx fy fz]; %#ok<AGROW>

            if elemCol > 0 && elemCol <= numel(tokens)
                elemSymbol = tokens{elemCol};
                % Clean symbol (remove numbers/charges like "Fe2+" -> "Fe")
                elemSymbol = regexprep(elemSymbol, '[^a-zA-Z]', '');
                thisMass = getAtomicMass(elemSymbol);
            else
                thisMass = 1.0;
            end
            asymmetricMasses = [asymmetricMasses; thisMass]; %#ok<AGROW>

        % --- Data row inside symmetry operations loop ---
        elseif inSymLoop
            % Preserve quoted strings containing spaces (e.g., 'x, y, z')
            symOpStr = strtrim(line);
            if strncmp(symOpStr, '_', 1) || strcmpi(symOpStr, 'loop_')
                inSymLoop = false;
                continue;
            end
            % Extract string token, removing leading/trailing quotes
            symOpStr = regexprep(symOpStr, '^[''"]|[''"]$', '');
            if ~isempty(symOpStr)
                symOps{end+1} = symOpStr; %#ok<AGROW>
            end
        end
    end
    fclose(fid);

    if isempty(a) || isempty(b) || isempty(c) || ...
       isempty(alpha) || isempty(beta) || isempty(gamma)
        error('parseCIFFile:missingCell', ...
            'Could not find all six cell parameters in %s.', filepath);
    end
    if isempty(asymmetricBasis)
        error('parseCIFFile:noAtoms', ...
            'No _atom_site_fract_x/y/z loop data found in %s.', filepath);
    end

    % --- Default to identity operator if no symmetry loop was present ---
    if isempty(symOps)
        symOps = {'x, y, z'};
    end

    % --- Apply symmetry operators & eliminate duplicates (special positions) ---
    basisFrac = [];
    basisMasses = [];
    tol = 1e-3; % Fractional tolerance for duplicate detection

    for i = 1:size(asymmetricBasis, 1)
        basePos = asymmetricBasis(i, :);
        thisMass = asymmetricMasses(i);
        for s = 1:numel(symOps)
            newPos = applySymOp(symOps{s}, basePos);
            newPos = mod(newPos, 1.0); % Wrap inside unit cell [0, 1)

            if ~isDuplicateAtom(basisFrac, newPos, tol)
                basisFrac = [basisFrac; newPos]; %#ok<AGROW>
                basisMasses = [basisMasses; thisMass]; %#ok<AGROW>
            end
        end
    end

    % --- Convert unit cell parameters to Cartesian lattice vectors ---
    alphaR = deg2rad(alpha);
    betaR  = deg2rad(beta);
    gammaR = deg2rad(gamma);

    a1 = [a, 0, 0];
    a2 = [b*cos(gammaR), b*sin(gammaR), 0];

    cx = c * cos(betaR);
    cy = c * (cos(alphaR) - cos(betaR)*cos(gammaR)) / sin(gammaR);
    czSq = c^2 - cx^2 - cy^2;
    if czSq < 0
        error('parseCIFFile:badAngles', ...
            'Cell angles in %s do not describe a valid unit cell.', filepath);
    end
    cz = sqrt(czSq);
    a3 = [cx, cy, cz];

    latticeVectors = [a1; a2; a3];

    fprintf('parseCIFFile: generated %d full cell basis atoms from %d asymmetric atoms in %s\n', ...
        size(basisFrac, 1), size(asymmetricBasis, 1), filepath);
end

% --- HELPER FUNCTIONS ---

function newPos = applySymOp(opStr, pos)
    % Evaluates a symmetry operation string safely
    x = pos(1); y = pos(2); z = pos(3); %#ok<NASGU>
    
    % Clean operation string and force lowercase for X, Y, Z
    opStr = lower(opStr);
    opStr = strrep(opStr, '''', '');
    opStr = strrep(opStr, '"', '');
    
    % Insert missing multiplication operators (e.g., '1/2x' becomes '1/2*x')
    opStr = regexprep(opStr, '(\d)\s*([xyz])', '$1*$2');
    
    parts = strsplit(opStr, ',');

    if numel(parts) ~= 3
        newPos = pos;
        return;
    end

    % Evaluate fractional arithmetic inside a try/catch block
    try
        newX = eval(parts{1});
        newY = eval(parts{2});
        newZ = eval(parts{3});
    catch ME
        warning('parseCIFFile:evalFailed', 'Could not evaluate symmetry operation "%s". Error: %s', opStr, ME.message);
        newX = x; newY = y; newZ = z;
    end

    newPos = [newX, newY, newZ];
end

function isDup = isDuplicateAtom(existingAtoms, newPos, tol)
    % Checks if newPos is identical to any position in existingAtoms,
    % accounting for unit cell boundary wraparound (e.g., 0.001 approx 0.999)
    if isempty(existingAtoms)
        isDup = false;
        return;
    end

    diffs = abs(bsxfun(@minus, existingAtoms, newPos));
    diffs = min(diffs, 1.0 - diffs); % Periodic boundary check
    dists = sqrt(sum(diffs.^2, 2));

    isDup = any(dists < tol);
end

function val = readCifNumber(line)
    tokens = strsplit(strtrim(line));
    tokens = tokens(~cellfun(@isempty, tokens));
    if numel(tokens) < 2
        error('parseCIFFile:badCellLine', 'Could not parse value from line: "%s"', line);
    end
    val = str2double(stripUncertainty(tokens{2}));
    if isnan(val)
        error('parseCIFFile:badCellLine', 'Could not parse a number from line: "%s"', line);
    end
end

function s = stripUncertainty(s)
    parenIdx = strfind(s, '(');
    if ~isempty(parenIdx)
        s = s(1:parenIdx(1)-1);
    end
end

function mass = getAtomicMass(symbol)
    persistent massMap
    if isempty(massMap)
        symbols = {'H','He','Li','Be','B','C','N','O','F','Ne', ...
             'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca', ...
             'Ti','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge', ...
             'As','Se','Br','Kr','Rb','Sr','Zr','Ag','Cd','In', ...
             'Sn','Sb','Te','I','Xe','Cs','Ba','Au','Hg','Pb','Bi'};
        massValues = [1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999, 18.998, 20.180, ...
             22.990, 24.305, 26.982, 28.085, 30.974, 32.06, 35.45, 39.948, 39.098, 40.078, ...
             47.867, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723, 72.630, ...
             74.922, 78.971, 79.904, 83.798, 85.468, 87.62, 91.224, 107.87, 112.41, 114.82, ...
             118.71, 121.76, 127.60, 126.90, 131.29, 132.91, 137.33, 196.97, 200.59, 207.2, 208.98];
        massMap = containers.Map(symbols, num2cell(massValues));
    end

    symbol = char(symbol);
    if isKey(massMap, symbol)
        mass = massMap(symbol);
    else
        warning('Element %s not in dictionary. Defaulting to 1.0 amu.', symbol);
        mass = 1.0;
    end
end