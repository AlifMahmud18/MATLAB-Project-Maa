function [candI, candJ, candDist, candU, shellDistances] = findPeriodicNeighborCandidates(pos, superLattice, imgRange)
%FINDPERIODICNEIGHBORCANDIDATES Every periodic-image neighbour pair (i<j) within the search range, vectorised.
%  [candI,candJ,candDist,candU,shellDistances] = findPeriodicNeighborCandidates(pos, superLattice, imgRange)
%  Returns pair indices, distances, unit vectors and the sorted distinct distances (= neighbour shells).
%  Self-images (i==j) are skipped: they add zero stiffness at the Gamma point.
%  Used by generateLatticeGeneral and autoRigidShells to pick bonds by shell.

    if nargin < 3 || isempty(imgRange)
        imgRange = -2:2;
    end

    Natoms = size(pos, 1);
    [sx, sy, sz] = ndgrid(imgRange, imgRange, imgRange);
    cartShifts = [sx(:), sy(:), sz(:)] * superLattice;
    nShifts = size(cartShifts, 1);

    [iiPairs, jjPairs] = find(triu(true(Natoms), 1));

    distCells = cell(nShifts, 1);
    uCells = cell(nShifts, 1);

    for s = 1:nShifts
        P2 = pos + cartShifts(s, :);

        dx = P2(jjPairs, 1) - pos(iiPairs, 1);
        dy = P2(jjPairs, 2) - pos(iiPairs, 2);
        dz = P2(jjPairs, 3) - pos(iiPairs, 3);
        d = sqrt(dx.^2 + dy.^2 + dz.^2);

        distCells{s} = d;
        uCells{s} = [dx ./ d, dy ./ d, dz ./ d];
    end

    candDist = vertcat(distCells{:});
    candU = vertcat(uCells{:});
    candI = repmat(iiPairs, nShifts, 1);
    candJ = repmat(jjPairs, nShifts, 1);

    roundedD = round(candDist * 1e4) / 1e4;
    shellDistances = sort(unique(roundedD));
end
