function [x, iters, converged] = gaussSeidelSolve(A, b, x0, tol, maxIter)
%GAUSSSEIDELSOLVE Iteratively solve A*x = b with the Gauss-Seidel method
%   (EEE 212 Exp. 5), always using the latest available values of x in
%   each update - exactly the "main point of the Gauss-Seidel iterative
%   process" called out in the labsheet.
%
%   [x, iters, converged] = gaussSeidelSolve(A, b, x0, tol, maxIter)
%
%   Convergence is guaranteed if A is strictly diagonally dominant, but
%   it also converges for any symmetric positive-definite A (the case
%   relevant to a spring-network stiffness matrix), even if A is not
%   diagonally dominant.
%
%   INPUTS
%     A       - [n x n] coefficient matrix
%     b       - [n x 1] right-hand side
%     x0      - (optional) [n x 1] initial guess. Default: zeros(n,1),
%               matching the labsheet's "very first iteration ... set
%               equal to zero" convention.
%     tol     - (optional) convergence tolerance on max|x_new - x_old|.
%               Default 1e-8.
%     maxIter - (optional) safety cap on iterations. Default 1000.
%
%   OUTPUTS
%     x         - [n x 1] approximate solution
%     iters     - number of iterations actually performed
%     converged - true if tol was reached before maxIter

    n = size(A, 1);
    if nargin < 3 || isempty(x0)
        x0 = zeros(n, 1);
    end
    if nargin < 4 || isempty(tol)
        tol = 1e-8;
    end
    if nargin < 5 || isempty(maxIter)
        maxIter = 1000;
    end

    x = x0(:);
    converged = false;

    for iters = 1:maxIter
        xOld = x;
        for i = 1:n
            s = A(i, :) * x - A(i, i) * x(i);   % uses latest x values throughout
            x(i) = (b(i) - s) / A(i, i);
        end
        if max(abs(x - xOld)) < tol
            converged = true;
            break;
        end
    end
end
