function [coords_sheared, force_distribution, Dmat_sheared, freqs_sheared, modes_sheared] = applyShearForce(coords, bonds, atomicMasses, K_matrix, springConstant, shearVal)
    %APPLYSHEARFORCE Shear the lattice (u_x = gamma*(z - z_min)) and re-solve its vibrations.
    %  shearVal is the engineering strain gamma itself, not a raw displacement.
    %  Colour = bond-force load per atom, sum of |k*(d - d0)| over its bonds (wrap-safe, from computeShearForces), scaled by the load at gamma = 0.5.
    %  Zero strain = zero colour. The net force is not used: it is exactly zero on centrosymmetric sites (MgO), so it would never change colour.
    %  Each bond unit vector is rotated by the affine map (wrap-around safe), then K and Dmat are rebuilt.
    %  Eigenvalues come from jacobiEigenSolver, sorted ascending; omega = sqrt(lambda) with no 2*pi.

    z_min = min(coords(:,3));

    displacement_field = zeros(size(coords));
    displacement_field(:,1) = shearVal * (coords(:,3) - z_min);
    coords_sheared = coords + displacement_field;

    S_xz = [0 0 1; 0 0 0; 0 0 0];  nAt = size(coords, 1);
    [~, ~, ~, fMag] = computeShearForces(coords, bonds, springConstant, shearVal, S_xz);
    [~, ~, ~, fRef] = computeShearForces(coords, bonds, springConstant, 0.5, S_xz);
    force_distribution = bondLoad(bonds, fMag, nAt);
    refMax = max(bondLoad(bonds, fRef, nAt));
    if refMax > 1e-12
        force_distribution = force_distribution / refMax;
    else
        force_distribution(:) = 0;
    end

    bonds_sheared = bonds;
    gamma_eff = shearVal;
    r_old = bonds(:, 3) .* bonds(:, 4:6);
    newVec = r_old;
    newVec(:, 1) = newVec(:, 1) + gamma_eff * r_old(:, 3);
    newLen = sqrt(sum(newVec.^2, 2));
    valid = newLen > 0;
    bonds_sheared(valid, 4:6) = newVec(valid, :) ./ newLen(valid);

    [K_sheared, M_matrix] = buildDynamicalMatrix(coords_sheared, bonds_sheared, springConstant, atomicMasses);

    invSqrtM = diag(1 ./ sqrt(diag(M_matrix)));
    Dmat_sheared = invSqrtM * K_sheared * invSqrtM;

    [modes_sheared, eigenvalues] = jacobiEigenSolver(Dmat_sheared);
    [eigenvalues, order] = sort(eigenvalues);
    modes_sheared = modes_sheared(:, order);
    freqs_sheared = sqrt(max(eigenvalues, 0));
end

function load = bondLoad(bonds, fMag, nAt)
    load = accumarray(bonds(:, 1), abs(fMag), [nAt 1]) + accumarray(bonds(:, 2), abs(fMag), [nAt 1]);
end
