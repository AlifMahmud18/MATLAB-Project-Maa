function [ratio, lamMin] = kSpaceRigidity(LV, frac, masses, cutoff, springConstant, nGrid)
%KSPACERIGIDITY Is the spring network stiff for waves at every k, not just at Gamma? (~1 healthy, ~0 floppy)
%   ratio = min/max over the Brillouin zone of lambda1(k)/|k|^2 (lambda1 = lowest eigenvalue of D(k)).
%   D(k) = M^-1/2 K(k) M^-1/2 from cellNeighborBonds; k on an nGrid^3 grid (default 6), folded into the first zone, Gamma excluded.
%   autoRigidShells accepts a shell count only if the Gamma rank test AND ratio >= 0.1 both pass.
%   Needed because the Gamma test alone can accept a network with near-zero modes at k ~= 0 (Al2O3 with 4 shells).

    if nargin < 6 || isempty(nGrid), nGrid = 6; end
    bonds = cellNeighborBonds(LV, frac, cutoff);
    nb = size(frac, 1);  dim = 3 * nb;  nB = size(bonds, 1);
    r = bonds(:, 3) .* bonds(:, 4:6);  u = bonds(:, 4:6);
    B = zeros(nB, 9);
    for a = 1:3, for b = 1:3
        B(:, 3*(a-1)+b) = springConstant * u(:, a) .* u(:, b);
    end, end
    bi = bonds(:, 1);  bj = bonds(:, 2);  [IA, IB, JA, JB] = deal(zeros(nB, 9));
    for a = 1:3, for b = 1:3
        c = 3 * (a - 1) + b;  IA(:, c) = 3*(bi-1) + a;  IB(:, c) = 3*(bi-1) + b;
        JA(:, c) = 3*(bj-1) + a;  JB(:, c) = 3*(bj-1) + b;
    end, end
    rows = [IA(:); JA(:); IA(:); JA(:)];  cols = [IB(:); JB(:); JB(:); IB(:)];
    Minv = kron(1 ./ sqrt(masses(:)), [1; 1; 1]);

    Brec = 2 * pi * inv(LV)';
    g = (-floor(nGrid/2):ceil(nGrid/2)-1) / nGrid;  [f1, f2, f3] = ndgrid(g);  F = [f1(:) f2(:) f3(:)];
    F(all(abs(F) < 1e-12, 2), :) = [];
    [d1, d2, d3] = ndgrid(-1:1);  Dl = [d1(:) d2(:) d3(:)];
    lam1 = zeros(size(F, 1), 1);  kk = zeros(size(F, 1), 1);
    for q = 1:size(F, 1)
        cand = (F(q, :) + Dl) * Brec;  [kk(q), ic] = min(vecnorm(cand, 2, 2));  kvec = cand(ic, :);
        ph = exp(1i * (r * kvec'));  Bp = B .* ph;  Bc = B .* conj(ph);
        H = full(sparse(rows, cols, [B(:); B(:); -Bp(:); -Bc(:)], dim, dim));
        ev = eig((Minv .* H .* Minv' + (Minv .* H .* Minv')') / 2, 'vector');
        lam1(q) = real(ev(1));
    end
    c = max(lam1, 0) ./ kk.^2;  lamMin = min(lam1);
    scale = springConstant / mean(masses) * min(bonds(:, 3))^2;
    if max(c) < 1e-4 * scale, ratio = 0; else, ratio = min(c) / max(c); end
end
