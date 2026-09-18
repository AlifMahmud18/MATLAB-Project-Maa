function [pos, bonds, masses, shellsUsed, numZero] = autoRigidShells(latVec, basis, basisMasses, N, springConstant, maxShells, requireKSpace)
%AUTORIGIDSHELLS Add neighbour shells one at a time until the spring network is rigid at Gamma AND at every k.
%  [pos,bonds,masses,shellsUsed,numZero] = autoRigidShells(latVec,basis,basisMasses,N,k,maxShells,requireKSpace)
%  Gamma test: rank(Dmat) leaves only the 3 translations. k-space test: kSpaceRigidity >= 0.1 (requireKSpace, default true).
%  Gamma alone can accept networks that are floppy for k ~= 0 (Al2O3: 4 shells fail, 6 pass; MgO 2, CaTiO3 10).
%  If no shell passes the k-space test the first Gamma-rigid shell is used; if maxShells (default 15) is hit, numZero > 3 shows it.

    if nargin < 6 || isempty(maxShells)
        maxShells = 15;
    end
    if nargin < 7 || isempty(requireKSpace)
        requireKSpace = true;
    end

    [pos, ~, masses] = generateLatticeGeneral(latVec, basis, basisMasses, N, 1);
    superLattice = N * latVec;

    [candI, candJ, candDist, candU, shellDistances] = ...
        findPeriodicNeighborCandidates(pos, superLattice);

    if isempty(shellDistances)
        error('autoRigidShells:noNeighbors', 'No neighbor interactions found in lattice.');
    end
    maxShells = min(maxShells, length(shellDistances));

    Natoms = size(pos, 1);
    totalDOF = 3 * Natoms;

    firstRigid = 0;
    for numShells = 1:maxShells
        cutoff = shellDistances(numShells) * 1.001;
        withinCutoff = candDist <= cutoff;
        bonds = [candI(withinCutoff), candJ(withinCutoff), candDist(withinCutoff), candU(withinCutoff, :)];

        [~, ~, Dmat] = buildDynamicalMatrix(pos, bonds, springConstant, masses);
        numZero = totalDOF - rank(Dmat, 1e-8);

        if numZero <= 3
            if firstRigid == 0, firstRigid = numShells; end
            % ASSUMED threshold 0.1 on the k-space rigidity ratio (healthy ~0.15-0.75, floppy ~0.001): a heuristic, not calibrated to a material.
            if ~requireKSpace || kSpaceRigidity(latVec, basis, basisMasses, cutoff, springConstant) >= 0.1
                shellsUsed = numShells;
                return;
            end
        end
    end

    if firstRigid > 0
        cutoff = shellDistances(firstRigid) * 1.001;  withinCutoff = candDist <= cutoff;
        bonds = [candI(withinCutoff), candJ(withinCutoff), candDist(withinCutoff), candU(withinCutoff, :)];
        numZero = 3;  shellsUsed = firstRigid;
        return;
    end

    shellsUsed = maxShells;
end
