function [pos, bonds, masses] = generateLattice(structureType, N, a, numShells)
%GENERATELATTICE Convenience wrapper around generateLatticeGeneral for the
%   built-in cubic structures ('sc', 'bcc', 'fcc').
%
%   [pos, bonds, masses] = generateLattice(structureType, N, a, numShells)
%
%   This is the function MainProject.m calls. It was missing from the
%   uploaded files (only generateLatticeGeneral.m, which takes explicit
%   lattice vectors/basis/masses, was present) - without it MainProject.m
%   errors with "Undefined function 'generateLattice'" on ANY MATLAB
%   version, not just R2016a. This wrapper just looks up the built-in
%   definition and forwards to generateLatticeGeneral, the same way
%   CrystalVibrationApp.m's BuildButtonPushed does.
%
%   NOTE ON 'a': as explained in getBuiltinCrystalDef.m, this project's
%   spring model only depends on bond DIRECTION, not length, so the
%   lattice constant 'a' never affects the computed frequencies - it's
%   accepted here only so MainProject.m's existing call signature keeps
%   working.

    if nargin < 3
        a = 1.0; %#ok<NASGU> % kept for interface compatibility only; see note above
    end
    if nargin < 4 || isempty(numShells)
        numShells = 1;
    end

    [latticeVectors, basisFrac, basisMasses] = getBuiltinCrystalDef(structureType);
    [pos, bonds, masses] = generateLatticeGeneral(latticeVectors, basisFrac, basisMasses, N, numShells);
end
