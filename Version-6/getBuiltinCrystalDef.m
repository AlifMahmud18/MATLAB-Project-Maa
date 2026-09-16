function [latticeVectors, basisFrac, basisMasses] = getBuiltinCrystalDef(name)
%GETBUILTINCRYSTALDEF Lattice vectors + basis for a few standard cubic
%   crystal structures, used as the GUI's default when no file is loaded.
%
%   [latticeVectors, basisFrac, basisMasses] = getBuiltinCrystalDef(name)
%
%   name (case-insensitive):
%     'sc'       - simple cubic, 1 atom/cell
%     'bcc'      - body-centered cubic, 2 atoms/cell
%     'fcc'      - face-centered cubic, 4 atoms/cell
%     'diamond'  - diamond cubic, 8 atoms/cell (two interpenetrating FCC
%                  sublattices, tetrahedral coordination - the classic
%                  example of a structure that is NOT locally
%                  centrosymmetric, so an imposed shear needs the
%                  internal-relaxation correction from
%                  relaxShearLattice.m)
%
%   All structures use a unit lattice constant a0 = 1 and a generic
%   atomic mass of 1.0 amu (scale as needed once loaded into the app).

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
