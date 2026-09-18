function [vonMisesYield, trescaYield, trescaStress] = ...
    checkPlasticity(principalStresses, vonMises, sigmaYield, tauYield)
%CHECKPLASTICITY Von Mises and Tresca yield tests for every atom's stress state.
%  [vonMisesYield, trescaYield, trescaStress] = checkPlasticity(principalStresses, vonMises, sigmaYield, tauYield)
%  Von Mises: yields if sigma_vm >= sigmaYield (distortion energy, ignores hydrostatic stress).
%  Tresca: yields if max(|s1-s2|,|s2-s3|,|s3-s1|)/2 >= tauYield (max shear, more conservative).
%  Fully vectorised; trescaStress is returned so callers can plot it.

    s1 = principalStresses(:,1);
    s2 = principalStresses(:,2);
    s3 = principalStresses(:,3);

    trescaStress = max([abs(s1 - s2), abs(s2 - s3), abs(s3 - s1)], [], 2) / 2;

    vonMisesYield = vonMises     >= sigmaYield;
    trescaYield   = trescaStress >= tauYield;
end
