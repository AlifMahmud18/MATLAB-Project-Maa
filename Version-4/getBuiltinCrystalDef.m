function [latticeVectors, basisFrac,basisMasses] = getBuiltinCrystalDef(structureType)
%GETBUILTINCRYSTALDEF Lattice vectors + fractional basis for the three
%   built-in cubic structures (sc, bcc, fcc), expressed as multiples of
%   a cubic lattice constant a = 1.
%
%   IMPORTANT PHYSICS NOTE: in this project's spring model, the restoring
%   force only depends on BOND DIRECTION (see buildDynamicalMatrix.m's
%   k*(u*u') formula) - it does NOT depend on bond length. That means
%   the lattice constant 'a' only changes the visual SIZE of the
%   structure, never the vibration frequencies. So there's no 'a'
%   parameter to tune here - it wouldn't do anything physically.

    switch lower(structureType)
        case 'sc'
            latticeVectors = eye(3);
            basisFrac = [0, 0, 0];
            basisMasses = 1.0; % Default 1.0 amu or element mass
        case 'bcc'
            latticeVectors = eye(3);
            basisFrac = [0,   0,   0;
                         0.5, 0.5, 0.5];
            basisMasses = [1.0; 1.0];

        case 'fcc'
            latticeVectors = eye(3);
            basisFrac = [0,   0,   0;
                         0.5, 0.5, 0;
                         0.5, 0,   0.5;
                         0,   0.5, 0.5];
            basisMasses = [1.0; 1.0; 1.0; 1.0];
        otherwise
            error('getBuiltinCrystalDef:badType', ...
                'structureType must be ''sc'', ''bcc'', or ''fcc''.');
    end
end
