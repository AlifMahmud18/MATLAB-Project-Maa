function [V, eigenvalues, sweeps] = jacobiEigenSolver(A, tol, maxSweeps)
%JACOBIEIGENSOLVER Compute eigenvalues/eigenvectors of a real symmetric
%   matrix using the classical CYCLIC JACOBI method. No built-in eig().
%
%   [V, eigenvalues, sweeps] = jacobiEigenSolver(A, tol, maxSweeps)
%
%   HOW THE JACOBI METHOD WORKS (the intuition)
%   A symmetric matrix A can be diagonalized by a sequence of simple 2D
%   rotations. Each rotation targets ONE off-diagonal pair (p,q) and
%   picks a rotation angle theta that makes A(p,q) exactly zero after
%   the rotation is applied (rotating both the corresponding row AND
%   column, since A must stay symmetric).
%
%   Doing this once for pair (p,q) can "un-zero" a different pair you
%   fixed earlier - but each full sweep (looping over every pair once)
%   provably shrinks the sum of squared off-diagonal elements. Repeating
%   sweeps drives the matrix to (near) diagonal form. The accumulated
%   product of all rotation matrices converges to the eigenvector matrix
%   V; the final diagonal entries are the eigenvalues.
%
%   THE ROTATION FORMULA (derived by requiring the new A(p,q) = 0):
%       theta = (A(q,q) - A(p,p)) / (2*A(p,q))
%       t     = sign(theta) / (abs(theta) + sqrt(theta^2 + 1))   [tan of rotation angle]
%       c     = 1 / sqrt(t^2 + 1)                                 [cos]
%       s     = t * c                                             [sin]
%
%   INPUTS
%     A         - [n x n] real symmetric matrix
%     tol       - (optional) convergence tolerance on off-diagonal norm.
%                 Default 1e-10.
%     maxSweeps - (optional) safety cap on number of sweeps. Default 100.
%
%   OUTPUTS
%     V           - [n x n] orthogonal matrix of eigenvectors (columns)
%     eigenvalues - [n x 1] eigenvalues, UNSORTED (matches column order of V)
%     sweeps      - number of sweeps actually performed (for diagnostics)
%
%   NOTE: eigenvalues come out in whatever order the algorithm leaves
%   them - sort them yourself afterward if you need ascending order:
%       [eigenvalues, order] = sort(eigenvalues);
%       V = V(:, order);

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
        % --- convergence check: sum of squares of off-diagonal elements ---
        offDiagNormSq = sum(A(:).^2) - sum(diag(A).^2);
        if sqrt(offDiagNormSq) < tol
            break;
        end

        % --- one full sweep: visit every (p,q) pair with p < q ---
        for p = 1:n-1
            for q = p+1:n
                apq = A(p, q);
                if abs(apq) < 1e-15
                    continue;   % already effectively zero, skip
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

                % --- update the 2x2 block (p,p),(p,q),(q,q) ---
                A(p, p) = app - t * apq;
                A(q, q) = aqq + t * apq;
                A(p, q) = 0;
                A(q, p) = 0;

                % --- update all other rows/columns i (i != p,q) ---
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

                % --- accumulate the rotation into V ---
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