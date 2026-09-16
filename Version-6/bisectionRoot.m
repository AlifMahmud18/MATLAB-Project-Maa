function [xMid, yMid, iters] = bisectionRoot(f, xLower, xUpper, xTol)
%BISECTIONROOT Bisection (interval-halving) root finder (EEE 212 Exp. 8).
%
%   [xMid, yMid, iters] = bisectionRoot(f, xLower, xUpper, xTol)
%
%   Implements the labsheet's "Algorithm: Bisection Method" verbatim:
%   input xLower, xUpper must bracket a root (f(xLower) and f(xUpper)
%   have opposite signs); the bracket is halved every iteration until
%   its width is <= xTol.
%
%   INPUTS
%     f      - function handle, f(x)
%     xLower - lower bracket
%     xUpper - upper bracket
%     xTol   - (optional) x-tolerance. Default 1e-8.
%
%   OUTPUTS
%     xMid  - approximate root
%     yMid  - f(xMid)
%     iters - number of bisections performed

    if nargin < 4 || isempty(xTol)
        xTol = 1e-8;
    end

    yLower = f(xLower);
    yUpper = f(xUpper);
    if yLower * yUpper > 0
        error('bisectionRoot:noBracket', ...
            'f(xLower) and f(xUpper) must have opposite signs.');
    end

    xMid = (xLower + xUpper) / 2.0;
    yMid = f(xMid);
    iters = 0;

    while (xUpper - xLower) / 2.0 > xTol
        iters = iters + 1;
        if yLower * yMid > 0
            xLower = xMid;
            yLower = yMid;
        else
            xUpper = xMid;
        end
        xMid = (xLower + xUpper) / 2.0;
        yMid = f(xMid);
    end
end
