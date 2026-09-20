function [stressTensor, principalStresses, vonMises] = ...
    computeAtomicStress(bonds, bondVecs, bondForceVec, Natoms, atomicVolume)
%COMPUTEATOMICSTRESS Per-atom virial (Clausius) stress tensor, principal
%   stresses, and von Mises stress for a central-force spring lattice.
%
%   [stressTensor, principalStresses, vonMises] = ...
%       computeAtomicStress(bonds, bondVecs, bondForceVec, Natoms, atomicVolume)
%
%   PHYSICS BACKGROUND (Study Guide Sec. 3)
%   ----------------------------------------------------------------
%   There is no literal "area" at the atomic scale, so molecular
%   simulation assigns each atom a well-defined stress using the VIRIAL
%   (Clausius / Tsai 1979) formula, built only from pairwise bond forces
%   and bond vectors:
%
%       sigma_i = (1/Omega_i) * sum_{j in bonded(i)} 1/2 * (f_ij (x) r_ij)
%
%   where (x) is the outer (tensor) product, Omega_i is the volume
%   attributed to atom i, f_ij is the force atom i feels from its bond
%   to atom j, and r_ij is the (current, deformed) vector from i to j.
%   Summed/averaged over a whole cell this reduces EXACTLY to the
%   macroscopic Cauchy stress tensor of continuum elasticity - here it
%   is kept per-atom to give a spatially resolved stress field.
%
%   SIGN CONVENTION: bondForceVec (from computeShearForces.m) is the
%   force ON atom i, and it points along +r_ij (toward atom j) whenever
%   the bond is in TENSION (stretched, bondForceMag > 0). So
%   outer(f_ij, r_ij) has a POSITIVE diagonal for a stretched bond -
%   i.e. tensile stress is positive and compressive stress is negative,
%   the standard continuum-mechanics convention.
%
%   Each bond is visited from BOTH ends: atom i accumulates
%   outer(+f, +r_ij), and atom j accumulates outer(-f, -r_ij), which is
%   the identical POSITIVE-for-tension quantity, using Newton's-third-
%   law force -f together with the opposite bond vector r_ji = -r_ij.
%   That symmetry is exactly why the code below adds the SAME "contrib"
%   term onto both endpoints instead of subtracting it for atom j.
%
%   PRINCIPAL STRESSES AND VON MISES
%   Diagonalizing each atom's 3x3 symmetric stress tensor gives the
%   principal stresses sigma1 >= sigma2 >= sigma3 - the normal stresses
%   along the axes where all shear components vanish. As the Study
%   Guide (Sec. 3.3) points out, this is mathematically IDENTICAL to
%   the 3x3 symmetric eigenproblem already solved for vibrational
%   modes, so jacobiEigenSolver.m (no built-in eig()) is reused
%   verbatim here, once per atom.
%
%   INPUTS
%     bonds        - [Nbonds x >=2] bond list; only columns 1 (i) and 2
%                     (j) are used
%     bondVecs     - [Nbonds x 3] CURRENT bond vectors r_j - r_i (output
%                     of computeShearForces.m)
%     bondForceVec - [Nbonds x 3] force ON atom i from each bond, +ve
%                     tension convention (output of computeShearForces.m)
%     Natoms       - number of atoms
%     atomicVolume - scalar OR [Natoms x 1] volume Omega_i attributed to
%                     each atom (e.g. unit-cell volume / atoms-per-cell
%                     for a monatomic lattice - see runStressStrainSweep.m)
%
%   OUTPUTS
%     stressTensor      - [3 x 3 x Natoms] full symmetric virial stress
%                          tensor for every atom
%     principalStresses - [Natoms x 3] principal stresses sigma1 >=
%                          sigma2 >= sigma3, one row per atom
%     vonMises           - [Natoms x 1] von Mises (distortion-energy)
%                          equivalent stress for every atom

    i_idx = bonds(:, 1);
    j_idx = bonds(:, 2);

    if isscalar(atomicVolume)
        Omega = atomicVolume * ones(Natoms, 1);
    else
        Omega = atomicVolume(:);
    end

    % ---- 1. Assemble the 6 independent tensor components, vectorized ----
    % Every bond contributes the SAME-sign 1/2*outer(f,r) term to both
    % of its endpoints (see the sign-convention note above), so this is
    % a single accumulation over the bond list - no explicit loop over
    % bonds or atoms is needed here.
    fx = bondForceVec(:,1); fy = bondForceVec(:,2); fz = bondForceVec(:,3);
    rx = bondVecs(:,1);     ry = bondVecs(:,2);     rz = bondVecs(:,3);

    % columns: xx, yy, zz, xy, yz, xz
    contrib = 0.5 * [fx.*rx, fy.*ry, fz.*rz, fx.*ry, fy.*rz, fx.*rz];

    S = zeros(Natoms, 6);
    for c = 1:6
        S(:, c) = accumarray(i_idx, contrib(:, c), [Natoms, 1]) ...
                + accumarray(j_idx, contrib(:, c), [Natoms, 1]);
    end
    S = S ./ Omega;   % divide by per-atom volume

    % ---- 2. Reassemble each atom's full symmetric 3x3 stress tensor ----
    stressTensor = zeros(3, 3, Natoms);
    stressTensor(1,1,:) = S(:,1); stressTensor(2,2,:) = S(:,2); stressTensor(3,3,:) = S(:,3);
    stressTensor(1,2,:) = S(:,4); stressTensor(2,1,:) = S(:,4);
    stressTensor(2,3,:) = S(:,5); stressTensor(3,2,:) = S(:,5);
    stressTensor(1,3,:) = S(:,6); stressTensor(3,1,:) = S(:,6);

    % ---- 3. Principal stresses via the required custom Jacobi solver ----
    % A true per-atom 3x3 eigendecomposition has no batched/vectorized
    % MATLAB equivalent (there is no "eig for a stack of matrices" that
    % avoids a loop), so this one loop - over atoms, calling
    % jacobiEigenSolver.m exactly as instructed, and never eig() - is a
    % deliberate, necessary exception to "avoid loops", not an oversight.
    principalStresses = zeros(Natoms, 3);
    for a = 1:Natoms
        [~, eigVals] = jacobiEigenSolver(stressTensor(:,:,a));
        principalStresses(a, :) = sort(eigVals, 'descend');
    end

    % ---- 4. Von Mises equivalent stress (Study Guide Sec. 3.3) ----
    s1 = principalStresses(:,1); s2 = principalStresses(:,2); s3 = principalStresses(:,3);
    vonMises = sqrt(0.5 * ((s1 - s2).^2 + (s2 - s3).^2 + (s3 - s1).^2));
end
