function dfdx = centralDiffDerivative(f, x, h)
%CENTRALDIFFDERIVATIVE Central difference approximation of f'(x), O(h^2)
%   (EEE 212 Exp. 6):
%
%       f'(x) = (f(x+h) - f(x-h)) / (2h) + O(h^2)
%
%   INPUTS
%     f - function handle, f(x)
%     x - point(s) at which to evaluate the derivative
%     h - (optional) step size. Default 1e-5.
%
%   OUTPUT
%     dfdx - approximate derivative, same size as x

    if nargin < 3 || isempty(h)
        h = 1e-5;
    end
    dfdx = (f(x + h) - f(x - h)) / (2 * h);
end
