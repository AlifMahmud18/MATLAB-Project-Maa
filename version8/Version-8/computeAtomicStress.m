function [stressTensor, principalStresses, vonMises] = ...
    computeAtomicStress(bonds, bondVecs, bondForceVec, Natoms, atomicVolume)
%COMPUTEATOMICSTRESS Per-atom virial stress tensor, principal stresses and von Mises stress.
%  [stressTensor, principalStresses, vonMises] = computeAtomicStress(bonds, bondVecs, bondForceVec, Natoms, atomicVolume)
%  sigma_i = (1/Omega_i) * sum over bonds of 1/2 * f (x) r, added on both ends of each bond; tension is positive.
%  Principal stresses come from jacobiEigenSolver on each 3x3 tensor (one loop over atoms).
%  von Mises = sqrt(1/2*((s1-s2)^2 + (s2-s3)^2 + (s3-s1)^2)).

    i_idx = bonds(:, 1);
    j_idx = bonds(:, 2);

    if isscalar(atomicVolume)
        Omega = atomicVolume * ones(Natoms, 1);
    else
        Omega = atomicVolume(:);
    end

    fx = bondForceVec(:,1); fy = bondForceVec(:,2); fz = bondForceVec(:,3);
    rx = bondVecs(:,1);     ry = bondVecs(:,2);     rz = bondVecs(:,3);

    contrib = 0.5 * [fx.*rx, fy.*ry, fz.*rz, fx.*ry, fy.*rz, fx.*rz];

    S = zeros(Natoms, 6);
    for c = 1:6
        S(:, c) = accumarray(i_idx, contrib(:, c), [Natoms, 1]) ...
                + accumarray(j_idx, contrib(:, c), [Natoms, 1]);
    end
    S = S ./ Omega;

    stressTensor = zeros(3, 3, Natoms);
    stressTensor(1,1,:) = S(:,1); stressTensor(2,2,:) = S(:,2); stressTensor(3,3,:) = S(:,3);
    stressTensor(1,2,:) = S(:,4); stressTensor(2,1,:) = S(:,4);
    stressTensor(2,3,:) = S(:,5); stressTensor(3,2,:) = S(:,5);
    stressTensor(1,3,:) = S(:,6); stressTensor(3,1,:) = S(:,6);

    principalStresses = zeros(Natoms, 3);
    for a = 1:Natoms
        [~, eigVals] = jacobiEigenSolver(stressTensor(:,:,a));
        principalStresses(a, :) = sort(eigVals, 'descend');
    end

    s1 = principalStresses(:,1); s2 = principalStresses(:,2); s3 = principalStresses(:,3);
    vonMises = sqrt(0.5 * ((s1 - s2).^2 + (s2 - s3).^2 + (s3 - s1).^2));
end
