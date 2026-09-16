function [root, iters] = newtonRaphsonRoot(f, x0, xTol, fDeriv)
%NEWTONRAPHSONROOT Newton-Raphson root finder (EEE 212 Exp. 8).
%
%   [root, iters] = newtonRaphsonRoot(f, x0, xTol, fDeriv)
%
%   Implements the labsheet's "Algorithm - Newton's Method" verbatim:
%       dx <- -f(x0) / f'(x0);  root <- x0 + dx
%       while abs(dx) > xTol
%           dx <- -f(root) / f'(root);  root <- root + dx
%       end
%
%   INPUTS
%     f      - function handle, f(x)
%     x0     - initial guess
%     xTol   - (optional) x-tolerance. Default 1e-8.
%     fDeriv - (optional) function handle for f'(x). If omitted, the
%              derivative is estimated with centralDiffDerivative.m
%              (Exp. 6), so this method also works when no analytic
%              derivative is available.
%
%   OUTPUTS
%     root  - approximate root
%     iters - number of Newton updates performed

    if nargin < 3 || isempty(xTol)
        xTol = 1e-8;
    end
    if nargin < 4 || isempty(fDeriv)
        fDeriv = @(x) centralDiffDerivative(f, x);
    end

    iters = 1;
    dx = -f(x0) / fDeriv(x0);
    root = x0 + dx;

    maxIter = 200;
    while abs(dx) > xTol && iters < maxIter
        dx = -f(root) / fDeriv(root);
        root = root + dx;
        iters = iters + 1;
    end
end
