function K = assembleShearHessian(bonds, r, kb, f, Natoms, includePrestress)
%ASSEMBLESHEARHESSIAN 3N x 3N Hessian of the bond network at the deformed bond vectors r.
%  Each bond adds the block  kb*u*u' + g*(I - u*u')  to (i,i),(j,j) and minus it to (i,j),(j,i).
%  kb = V'' (bond stiffness); g = f/d is the pre-stress term, used only if includePrestress.
%  Broken bonds (kb = f = 0) drop out. Built with sparse() so duplicate entries add.

    d = sqrt(sum(r.^2, 2));
    u = r ./ d;
    if includePrestress
        g = f ./ d;
    else
        g = zeros(size(d));
    end

    i = bonds(:, 1);  j = bonds(:, 2);
    nB = numel(i);
    I = zeros(36 * nB, 1);  J = I;  V = I;
    p = 0;
    for a = 1:3
        for b = 1:3
            Bab = (kb - g) .* u(:, a) .* u(:, b) + g * (a == b);
            ia = 3*(i-1) + a;  ib = 3*(i-1) + b;
            ja = 3*(j-1) + a;  jb = 3*(j-1) + b;
            rows = [ia; ja; ia; ja];
            cols = [ib; jb; jb; ib];
            vals = [Bab; Bab; -Bab; -Bab];
            n = numel(rows);
            I(p+1:p+n) = rows;  J(p+1:p+n) = cols;  V(p+1:p+n) = vals;
            p = p + n;
        end
    end
    K = full(sparse(I(1:p), J(1:p), V(1:p), 3*Natoms, 3*Natoms));
    K = (K + K') / 2;
end
