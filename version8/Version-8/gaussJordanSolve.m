function x = gaussJordanSolve(A, b)
%GAUSSJORDANSOLVE Solve A*x = b by Gaussian elimination with partial pivoting and back substitution (EEE 212 Exp. 5).
%  x = gaussJordanSolve(A, b)   A: n x n, b: n x m. Raises 'gaussJordanSolve:singular' if a pivot is ~0.
%  Used by relaxShearLattice (interior-atom equilibrium) and runStressStrainSweep (least-squares fit).

    n = size(A, 1);
    if size(A, 2) ~= n
        error('gaussJordanSolve:notSquare', 'A must be square.');
    end
    if size(b, 1) ~= n
        error('gaussJordanSolve:sizeMismatch', 'b must have %d rows.', n);
    end

    aug = [A, b];
    nRhs = size(b, 2);

    for k = 1:n-1
        [maxVal, mRel] = max(abs(aug(k:n, k)));
        p = mRel + k - 1;
        if maxVal < 1e-14
            error('gaussJordanSolve:singular', ...
                'Matrix is singular (or nearly so) at column %d.', k);
        end
        if p ~= k
            aug([k p], :) = aug([p k], :);
        end

        for i = k+1:n
            u = aug(i, k) / aug(k, k);
            aug(i, k:end) = aug(i, k:end) - u * aug(k, k:end);
        end
    end

    if abs(aug(n, n)) < 1e-14
        error('gaussJordanSolve:singular', 'Matrix is singular (or nearly so) at the last pivot.');
    end

    x = zeros(n, nRhs);
    x(n, :) = aug(n, n+1:end) / aug(n, n);
    for i = n-1:-1:1
        s = aug(i, i+1:n) * x(i+1:n, :);
        x(i, :) = (aug(i, n+1:end) - s) / aug(i, i);
    end
end
