function bonds = cellNeighborBonds(LV, frac, cutoff)
%CELLNEIGHBORBONDS Unit-cell bond list [i j d0 ux uy uz] within cutoff, all periodic images included.
%   Bonds between an atom and its own periodic image (i == j) are kept, each counted once.
%   They add nothing at Gamma but do matter for waves with k ~= 0, so generateLatticeGeneral's list is not enough for k-space work.
%   Used by kSpaceRigidity.m and brillouinShearSurface.m.

    nb = size(frac, 1);  P = frac * LV;
    nmax = min(4, ceil(cutoff / min(vecnorm(LV, 2, 2))) + 1);
    [n1, n2, n3] = ndgrid(-nmax:nmax);  Nn = [n1(:) n2(:) n3(:)];  T = Nn * LV;
    half = Nn(:, 1) > 0 | (Nn(:, 1) == 0 & (Nn(:, 2) > 0 | (Nn(:, 2) == 0 & Nn(:, 3) > 0)));
    bonds = zeros(0, 6);
    for i = 1:nb
        for j = i:nb
            R = P(j, :) - P(i, :) + T;  d = sqrt(sum(R.^2, 2));
            ok = d > 1e-9 & d <= cutoff * 1.0005;
            if i == j, ok = ok & half; end
            bonds = [bonds; repmat([i j], nnz(ok), 1), d(ok), R(ok, :) ./ d(ok)]; %#ok<AGROW>
        end
    end
end
