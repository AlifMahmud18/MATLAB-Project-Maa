function [vonMisesYield, trescaYield, trescaStress] = ...
    checkPlasticity(principalStresses, vonMises, sigmaYield, tauYield)
%CHECKPLASTICITY Evaluate the Von Mises and Tresca yield criteria for
%   every atom's stress state.
%
%   [vonMisesYield, trescaYield, trescaStress] = ...
%       checkPlasticity(principalStresses, vonMises, sigmaYield, tauYield)
%
%   PHYSICS BACKGROUND (Study Guide Sec. 4)
%   ----------------------------------------------------------------
%   Both criteria take the SAME input - the three principal stresses -
%   and convert a full 3D stress state into a single pass/fail number;
%   that is why this is one function with two outputs rather than two
%   separate functions, as the roadmap specifies.
%
%   VON MISES (distortion-energy, Sec. 4.2) - the default/first choice
%   for ductile, isotropic materials such as metals. Physically:
%   yielding occurs when the elastic distortion (shape-change) energy
%   density reaches a critical value; it does NOT depend on hydrostatic
%   (volume-changing) stress, matching the observed behavior of ductile
%   metals under pure hydrostatic pressure.
%
%       yields if  sigma_vm >= sigma_yield
%       sigma_vm = sqrt( 1/2 * [(s1-s2)^2 + (s2-s3)^2 + (s3-s1)^2] )
%
%   TRESCA (maximum shear stress, Sec. 4.3) - motivated by the physical
%   picture that crystalline slip is a shear-driven process
%   (dislocation glide). It is MORE CONSERVATIVE than von Mises (predicts
%   yield at a slightly lower load) for most stress states:
%
%       yields if  max(|s1-s2|, |s2-s3|, |s3-s1|) / 2 >= tau_yield
%
%   INPUTS
%     principalStresses - [Natoms x 3] sigma1 >= sigma2 >= sigma3 per
%                          atom (from computeAtomicStress.m)
%     vonMises           - [Natoms x 1] von Mises stress per atom (from
%                          computeAtomicStress.m)
%     sigmaYield         - scalar Von Mises yield strength sigma_yield
%     tauYield           - scalar Tresca (max shear) yield strength
%                          tau_yield
%
%   OUTPUTS
%     vonMisesYield - [Natoms x 1] logical, true where sigma_vm >= sigmaYield
%     trescaYield   - [Natoms x 1] logical, true where the Tresca
%                      maximum-shear stress >= tauYield
%     trescaStress  - [Natoms x 1] the Tresca maximum-shear-stress value
%                      itself, returned so callers (e.g. the
%                      stress-strain sweep) can plot/report it alongside
%                      von Mises without recomputing it

    s1 = principalStresses(:,1);
    s2 = principalStresses(:,2);
    s3 = principalStresses(:,3);

    % Fully vectorized: all three candidate shear values for every atom
    % at once, then the row-wise max - no loop over atoms needed.
    trescaStress = max([abs(s1 - s2), abs(s2 - s3), abs(s3 - s1)], [], 2) / 2;

    vonMisesYield = vonMises     >= sigmaYield;
    trescaYield   = trescaStress >= tauYield;
end
