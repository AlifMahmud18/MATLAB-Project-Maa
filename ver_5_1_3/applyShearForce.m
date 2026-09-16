function [coords_sheared, force_distribution, Dmat_sheared, freqs_sheared, modes_sheared] = applyShearForce(coords, bonds, atomicMasses, K_matrix, springConstant, shearVal)
    % APPLYSHEARFORCE Applies a vertical shear strain profile across the crystal unit cell.
    %
    %   buildDynamicalMatrix needs the bond list and the scalar spring
    %   constant to rebuild K for the sheared geometry - it has no
    %   "elementTypes" input - so those are what this function now takes
    %   in place of the old, unused elementTypes argument.
    
    % 1. Identify coordinate bounds for shear application
    y_min = min(coords(:,2));
    y_max = max(coords(:,2));
    span_y = max(y_max - y_min, eps);
    
    % 2. Compute displaced coordinates 
    coords_sheared = coords;
    displacement_field = zeros(size(coords));
    
    for i = 1:size(coords, 1)
        factor = (coords(i,2) - y_min) / span_y;
        dx = shearVal * factor; 
        coords_sheared(i,1) = coords(i,1) + dx;
        displacement_field(i,1) = dx;
    end
    
    % 3. Calculate local force distribution magnitude
    num_atoms = size(coords, 1);
    force_distribution = zeros(num_atoms, 1);
    for i = 1:num_atoms
        idx_range = (3*i-2):(3*i);
        u_local = reshape(displacement_field(i,:), [], 1);
        K_local = K_matrix(idx_range, idx_range);
        force_distribution(i) = norm(K_local * u_local);
    end
    
    % Normalize force distribution for colormap scaling [0, 1]
    if max(force_distribution) > 0
        force_distribution = force_distribution / max(force_distribution);
    end
    
    % 4. Rebuild dynamical matrix for the sheared configuration
    % Note: Make sure buildDynamicalMatrix is in your working path!
    % buildDynamicalMatrix's real signature is
    % (pos, bonds, springConstant, masses) - NOT (pos, elementTypes, masses).
    % IMPORTANT: buildDynamicalMatrix only uses "pos" to count atoms - it
    % builds every bond's stiffness block from the unit vector already
    % stored in bonds(:,4:6), and never recomputes that direction from
    % pos itself. So if we pass the ORIGINAL bonds list here, K_sheared
    % comes out identical for every shearVal and the vibrational modes
    % never move. Fix: re-derive each bond's unit vector from the
    % sheared coordinates before rebuilding K.
    bonds_sheared = bonds;
    for b = 1:size(bonds, 1)
        i = bonds(b, 1);
        j = bonds(b, 2);
        dvec = coords_sheared(j, :) - coords_sheared(i, :);
        len = norm(dvec);
        if len > 0
            bonds_sheared(b, 4:6) = dvec / len;
        end
    end

    [K_sheared, M_matrix] = buildDynamicalMatrix(coords_sheared, bonds_sheared, springConstant, atomicMasses);
    
    invSqrtM = diag(1 ./ sqrt(diag(M_matrix)));
    Dmat_sheared = invSqrtM * K_sheared * invSqrtM;
    
    % 5. Solve eigenproblem using your custom Jacobi solver
    % jacobiEigenSolver returns [V, eigenvalues, sweeps] - V first, then
    % eigenvalues - so grab them in that order (they were swapped here,
    % which meant "eigenvalues" below was actually the eigenvector matrix).
    %
    % jacobiEigenSolver's own docs say its eigenvalues come out UNSORTED,
    % while the main window's spectrum (app.omega) IS sorted ascending -
    % so without sorting here too, this plot's mode order won't match the
    % main spectrum's shape at all, even for the pristine (shearVal=0)
    % case. Also: app.omega = sqrt(eigVals) directly (angular frequency,
    % omega) with no /(2*pi) - dividing by 2*pi here would silently
    % rescale this plot relative to the main spectrum. Match both so the
    % two are directly comparable.
    [modes_sheared, eigenvalues] = jacobiEigenSolver(Dmat_sheared);
    [eigenvalues, order] = sort(eigenvalues);
    modes_sheared = modes_sheared(:, order);
    freqs_sheared = sqrt(max(eigenvalues, 0));
end