function dfdx = centralDiffDerivative(f, x, h)
%CENTRALDIFFDERIVATIVE f'(x) ~ (f(x+h) - f(x-h)) / (2h), error O(h^2) (EEE 212 Exp. 6).
%  dfdx = centralDiffDerivative(f, x, h)   f: function handle, x: point(s), h: step (default 1e-5).

    if nargin < 3 || isempty(h)
        h = 1e-5;
    end
    dfdx = (f(x + h) - f(x - h)) / (2 * h);
end
