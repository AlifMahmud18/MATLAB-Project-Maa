function [I, Icum] = compositeTrapz(y, h)
%COMPOSITETRAPZ Composite trapezoidal rule for tabulated data on a
%   uniform grid (EEE 212 Exp. 7).
%
%   [I, Icum] = compositeTrapz(y, h)
%
%   Matches the labsheet's algorithm for a tabulated function:
%       sum <- (f1 + f_(n+1)) / 2
%       sum <- sum + f_j        for j = 2 .. n
%       integral <- h * sum
%   computed here cumulatively (Icum) so the running integral is
%   available at every sample point, not just the final total;
%   Icum(end) equals exactly the labsheet's single "integral" result.
%
%   INPUTS
%     y - tabulated values f(x1), f(x1+h), ..., f(x1+n*h)
%     h - uniform step size
%
%   OUTPUTS
%     I    - total integral over the whole tabulated range
%     Icum - cumulative integral at each sample point (Icum(1) = 0)

    y = y(:);
    n = numel(y);
    Icum = zeros(n, 1);
    for k = 2:n
        Icum(k) = Icum(k-1) + h * (y(k) + y(k-1)) / 2;
    end
    I = Icum(end);
end
