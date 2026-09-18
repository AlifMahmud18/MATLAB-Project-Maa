function gap = shearYieldGap(gamma, pos, bonds, springConstant, Natoms, atomicVolume, criterion, threshold)
%SHEARYIELDGAP Root-finding target: worst-atom stress minus the yield threshold at strain gamma.
%  gap = shearYieldGap(gamma, pos, bonds, k, Natoms, atomicVolume, criterion, threshold)
%  criterion is 'vonmises' or 'tresca'; gap = 0 at the yield strain (bisection/Newton in runStressStrainSweep).
%  Re-runs computeShearForces and computeAtomicStress at the given gamma.

    [~, ~, bondVecs, ~, bondForceVec, ~] = computeShearForces(pos, bonds, springConstant, gamma);
    [~, principalStresses, vonMises] = computeAtomicStress(bonds, bondVecs, bondForceVec, Natoms, atomicVolume);

    switch lower(criterion)
        case 'vonmises'
            gap = max(vonMises) - threshold;
        case 'tresca'
            [~, ~, trescaStress] = checkPlasticity(principalStresses, vonMises, Inf, Inf);
            gap = max(trescaStress) - threshold;
        otherwise
            error('shearYieldGap:badCriterion', 'criterion must be ''vonmises'' or ''tresca''.');
    end
end
