function [coords_sheared, dispField, bondVecs_sheared, bondForceMag, bondForceVec, atomForce] = ...
    computeShearForces(pos, bonds, springConstant, gamma, S)
%COMPUTESHEARFORCES Affine (Cauchy-Born) shear u = gamma*S*r and the Hooke's-law bond forces it causes.
%  [coords_sheared,dispField,bondVecs_sheared,bondForceMag,bondForceVec,atomForce] = computeShearForces(pos,bonds,k,gamma,S)
%  Bond vectors are transformed directly, r_new = (I + gamma*S) r_old, so wrapped bonds stay correct.
%  f = k*(d_new - d0), positive = tension; the force on atom i points along the bond, on atom j the opposite.
%  S defaults to the xy shear [0 1 0; 0 0 0; 0 0 0]; atomForce is the net force on each atom.

    if nargin < 5 || isempty(S)
        S = [0 1 0; 0 0 0; 0 0 0];
    end

    Natoms = size(pos, 1);

    dispField      = gamma * (pos * S');
    coords_sheared = pos + dispField;

    i_idx = bonds(:, 1);
    j_idx = bonds(:, 2);
    d0    = bonds(:, 3);
    r_old = d0 .* bonds(:, 4:6);

    bondVecs_sheared = r_old + gamma * (r_old * S');
    d_new            = sqrt(sum(bondVecs_sheared.^2, 2));
    unit_new         = bondVecs_sheared ./ d_new;

    bondForceMag = springConstant * (d_new - d0);
    bondForceVec = bondForceMag .* unit_new;

    atomForce = zeros(Natoms, 3);
    for d = 1:3
        atomForce(:, d) = accumarray(i_idx, bondForceVec(:, d), [Natoms, 1]) ...
                         - accumarray(j_idx, bondForceVec(:, d), [Natoms, 1]);
    end
end
