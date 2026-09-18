function [latticeVectors, basisFrac, basisMasses] = getBuiltinCrystalDef(name)
%GETBUILTINCRYSTALDEF Lattice vectors, fractional basis and masses of 'sc', 'bcc', 'fcc' or 'diamond' (a0 = 1, mass = 1).
%  [latticeVectors, basisFrac, basisMasses] = getBuiltinCrystalDef(name)   name is case-insensitive.
%  Default crystal of the GUI when no file is loaded.

    % Reduced units (a0 = 1, mass = 1): scale to the real lattice constant and atomic mass for a real material.
    a0 = 1.0;
    mass = 1.0;

    switch lower(name)
        case 'sc'
            latticeVectors = a0 * eye(3);
            basisFrac = [0 0 0];
            basisMasses = mass;

        case 'bcc'
            latticeVectors = a0 * eye(3);
            basisFrac = [0 0 0; 0.5 0.5 0.5];
            basisMasses = [mass; mass];

        case 'fcc'
            latticeVectors = a0 * eye(3);
            basisFrac = [0 0 0; 0.5 0.5 0; 0.5 0 0.5; 0 0.5 0.5];
            basisMasses = mass * ones(4, 1);

        case 'diamond'
            latticeVectors = a0 * eye(3);
            basisFrac = [0 0 0; 0.5 0.5 0; 0.5 0 0.5; 0 0.5 0.5; ...
                         0.25 0.25 0.25; 0.75 0.75 0.25; 0.75 0.25 0.75; 0.25 0.75 0.75];
            basisMasses = mass * ones(8, 1);

        otherwise
            error('getBuiltinCrystalDef:unknownName', ...
                'Unknown builtin crystal "%s". Choose sc, bcc, fcc, or diamond.', name);
    end
end
