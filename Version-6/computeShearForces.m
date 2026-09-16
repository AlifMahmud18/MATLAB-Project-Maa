function [coords_sheared, dispField, bondVecs_sheared, bondForceMag, bondForceVec, atomForce] = ...
    computeShearForces(pos, bonds, springConstant, gamma, S)
%COMPUTESHEARFORCES Apply an affine (Cauchy-Born) shear deformation to the
%   lattice and compute the resulting bond-by-bond Hooke's-law restoring
%   forces and the net force on every atom.
%
%   [coords_sheared, dispField, bondVecs_sheared, bondForceMag, ...
%       bondForceVec, atomForce] = computeShearForces(pos, bonds, ...
%       springConstant, gamma, S)
%
%   PHYSICS BACKGROUND (Study Guide Sec. 2)
%   ----------------------------------------------------------------
%   An AFFINE (homogeneous) deformation applies the same linear map to
%   every atom regardless of its position - this is the Cauchy-Born
%   rule, the standard assumption connecting an imposed continuum strain
%   field to individual atomic displacements in a defect-free crystal:
%
%       u(r) = gamma * S * r          (r = equilibrium position vector)
%
%   For a pure engineering shear in the x-y plane, S has S(1,2) = 1 and
%   every other entry zero, which reduces to the familiar
%
%       u_x = gamma * y,   u_y = 0,   u_z = 0
%
%   Once every atom has moved to its sheared position, each bond's
%   length has changed. Because our spring model is CENTRAL-FORCE (the
%   same assumption already baked into buildDynamicalMatrix.m), the
%   restoring force on a stretched/compressed bond is plain 1D Hooke's
%   law applied along the CURRENT (deformed) bond axis:
%
%       f = k * (d_new - d_0)                  (scalar, +ve = tension)
%       f_vec_on_i = f * (r_j_new - r_i_new) / d_new
%
%   A positive f pulls atom i toward atom j (the bond wants to shrink
%   back to d_0), and by Newton's third law atom j feels exactly
%   -f_vec_on_i. Summing every bond's contribution at each atom (a
%   discrete divergence of the bond-force field) gives the net atomic
%   force. This is exactly zero in a perfectly homogeneous infinite
%   lattice at equilibrium, but is generally NON-zero near the free
%   surfaces/edges of a finite supercell - that residual is a real,
%   physically meaningful "surface relaxation" force that this
%   idealized affine model does not relax away (Study Guide Sec. 2.4,
%   "affine vs. relaxation").
%
%   INPUTS
%     pos            - [Natoms x 3] equilibrium atomic positions
%     bonds          - [Nbonds x 6] bond list (i, j, d0, ux0, uy0, uz0),
%                       exactly the format produced by
%                       generateLatticeGeneral.m / consumed by
%                       buildDynamicalMatrix.m. Only columns 1, 2, and 3
%                       (i, j, d0) are used here - the ORIGINAL unit
%                       vector in columns 4:6 is intentionally ignored,
%                       because after shearing the bond direction must
%                       be recomputed from the new, deformed positions
%                       (see applyShearForce.m for a worked example of
%                       exactly this pitfall, and why it matters).
%     springConstant - scalar spring constant k (same convention as
%                       buildDynamicalMatrix.m)
%     gamma          - scalar engineering shear strain to apply
%     S              - (optional) [3x3] shear direction tensor. Default
%                       is the standard xy engineering shear,
%                       S = [0 1 0; 0 0 0; 0 0 0], i.e. u_x = gamma*y.
%
%   OUTPUTS
%     coords_sheared    - [Natoms x 3] atom positions after the affine
%                          shear (r_old + gamma*S*r_old)
%     dispField         - [Natoms x 3] displacement u = coords_sheared - pos
%     bondVecs_sheared  - [Nbonds x 3] CURRENT (deformed) bond vectors,
%                          r_j_new - r_i_new, one row per bond - reused
%                          directly by computeAtomicStress.m as the r_ij
%                          in the virial formula
%     bondForceMag      - [Nbonds x 1] scalar bond force, +ve = tension
%                          (bond stretched, pulling the two atoms
%                          together), -ve = compression (pushing apart)
%     bondForceVec      - [Nbonds x 3] force vector ON ATOM i (the first
%                          index in that bond's row) from that bond;
%                          atom j feels exactly -bondForceVec(row,:).
%     atomForce         - [Natoms x 3] net Hooke's-law force on every
%                          atom, summed over all of its bonds

    if nargin < 5 || isempty(S)
        S = [0 1 0; 0 0 0; 0 0 0];   % default: pure xy engineering shear
    end

    Natoms = size(pos, 1);

    % ---- 1. Affine (Cauchy-Born) displacement, fully vectorized ----
    % u = gamma * S * r for every atom at once: (S*pos')' = pos*S'
    dispField      = gamma * (pos * S');
    coords_sheared = pos + dispField;

    % ---- 2. Recompute every bond's CURRENT vector and length ----
    % (must be re-derived from the sheared geometry - the pre-shear unit
    % vector in bonds(:,4:6) goes stale the moment gamma ~= 0, exactly
    % the bug flagged and fixed in applyShearForce.m)
    i_idx = bonds(:, 1);
    j_idx = bonds(:, 2);
    d0    = bonds(:, 3);

    bondVecs_sheared = coords_sheared(j_idx, :) - coords_sheared(i_idx, :);
    d_new            = sqrt(sum(bondVecs_sheared.^2, 2));
    unit_new         = bondVecs_sheared ./ d_new;       % [Nbonds x 3]

    % ---- 3. Bond-by-bond central-force Hooke's law (vectorized) ----
    bondForceMag = springConstant * (d_new - d0);       % [Nbonds x 1]
    bondForceVec = bondForceMag .* unit_new;            % force ON atom i

    % ---- 4. Discrete divergence: sum bond forces onto each atom ----
    % accumarray replaces what would otherwise be a "for each bond, add
    % to atom i, subtract from atom j" loop over Nbonds. The loop below
    % only runs 3 times (once per x/y/z column), never once per atom or
    % bond, so this stays fully vectorized in the sense that matters.
    % Atom j gets the reaction force -bondForceVec by Newton's third law.
    atomForce = zeros(Natoms, 3);
    for d = 1:3
        atomForce(:, d) = accumarray(i_idx, bondForceVec(:, d), [Natoms, 1]) ...
                         - accumarray(j_idx, bondForceVec(:, d), [Natoms, 1]);
    end
end
