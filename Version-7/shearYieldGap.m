function gap = shearYieldGap(gamma, pos, bonds, springConstant, Natoms, atomicVolume, criterion, threshold)
%SHEARYIELDGAP Scalar root-finding target for locating the exact yield
%   strain (EEE 212 Exp. 8): how far the worst-stressed atom is from a
%   yield threshold at a given engineering shear strain gamma.
%
%   gap = shearYieldGap(gamma, pos, bonds, springConstant, Natoms, ...
%             atomicVolume, criterion, threshold)
%
%   gap(gamma) = 0 exactly at the onset of yielding, so it is the
%   function handed to bisectionRoot.m / newtonRaphsonRoot.m to refine
%   the coarse sweep-sample-based yield estimate in runStressStrainSweep.m
%   into a precise strain value.
%
%   INPUTS
%     gamma          - scalar engineering shear strain to evaluate at
%     pos, bonds     - pristine lattice (from generateLatticeGeneral.m)
%     springConstant - scalar spring constant k
%     Natoms         - number of atoms
%     atomicVolume   - scalar or [Natoms x 1] per-atom volume
%     criterion      - 'vonmises' or 'tresca'
%     threshold      - sigmaYield (von Mises) or tauYield (Tresca)
%
%   OUTPUT
%     gap - max_i(stress_i) - threshold, using the same "does ANY atom
%           yield" criterion as checkPlasticity.m

    [~, ~, bondVecs, ~, bondForceVec, ~] = computeShearForces(pos, bonds, springConstant, gamma);
    [~, principalStresses, vonMises] = computeAtomicStress(bonds, bondVecs, bondForceVec, Natoms, atomicVolume);

    switch lower(criterion)
        case 'vonmises'
            gap = max(vonMises) - threshold;
        case 'tresca'
            % Inf thresholds: only the trescaStress output is wanted here.
            [~, ~, trescaStress] = checkPlasticity(principalStresses, vonMises, Inf, Inf);
            gap = max(trescaStress) - threshold;
        otherwise
            error('shearYieldGap:badCriterion', 'criterion must be ''vonmises'' or ''tresca''.');
    end
end
