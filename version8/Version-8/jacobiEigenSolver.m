function [V, eigenvalues, sweeps] = jacobiEigenSolver(A, tol, maxSweeps)
%JACOBIEIGENSOLVER Eigenvalues and eigenvectors of a real symmetric matrix by the cyclic Jacobi method (no eig()).
%  [V, eigenvalues, sweeps] = jacobiEigenSolver(A, tol, maxSweeps)   defaults tol = 1e-10, maxSweeps = 100.
%  Each rotation zeroes one off-diagonal pair (p,q); sweeps repeat until the off-diagonal norm < tol.
%  Eigenvalues are UNSORTED (same order as the columns of V); sort them if needed.

    if nargin < 2 || isempty(tol)
        tol = 1e-10;
    end
    if nargin < 3 || isempty(maxSweeps)
        maxSweeps = 100;
    end

    n = size(A, 1);
    if size(A, 2) ~= n
        error('jacobiEigenSolver:notSquare', 'A must be square.');
    end
    if norm(A - A', 'fro') > 1e-8 * norm(A, 'fro')
        warning('jacobiEigenSolver:notSymmetric', ...
            'A does not look symmetric - results may be wrong.');
    end

    V = eye(n);

    for sweeps = 1:maxSweeps
        offDiagNormSq = sum(A(:).^2) - sum(diag(A).^2);
        if sqrt(offDiagNormSq) < tol
            break;
        end

        for p = 1:n-1
            for q = p+1:n
                apq = A(p, q);
                if abs(apq) < 1e-15
                    continue;
                end

                app = A(p, p);
                aqq = A(q, q);

                theta = (aqq - app) / (2 * apq);
                if theta >= 0
                    t = 1 / (theta + sqrt(theta^2 + 1));
                else
                    t = -1 / (-theta + sqrt(theta^2 + 1));
                end
                c = 1 / sqrt(t^2 + 1);
                s = t * c;

                A(p, p) = app - t * apq;
                A(q, q) = aqq + t * apq;
                A(p, q) = 0;
                A(q, p) = 0;

                for i = 1:n
                    if i ~= p && i ~= q
                        aip = A(i, p);
                        aiq = A(i, q);
                        A(i, p) = c * aip - s * aiq;
                        A(p, i) = A(i, p);
                        A(i, q) = s * aip + c * aiq;
                        A(q, i) = A(i, q);
                    end
                end

                for i = 1:n
                    vip = V(i, p);
                    viq = V(i, q);
                    V(i, p) = c * vip - s * viq;
                    V(i, q) = s * vip + c * viq;
                end
            end
        end
    end

    eigenvalues = diag(A);

    if sweeps == maxSweeps
        warning('jacobiEigenSolver:notConverged', ...
            'Reached maxSweeps (%d) without hitting tolerance %.2e.', maxSweeps, tol);
    end
end
