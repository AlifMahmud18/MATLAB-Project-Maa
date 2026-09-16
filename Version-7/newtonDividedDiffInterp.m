function yq = newtonDividedDiffInterp(x, y, xq)
%NEWTONDIVIDEDDIFFINTERP Newton polynomial interpolation via divided
%   differences (EEE 212 Exp. 3), evaluated with nested (Horner-style)
%   multiplication.
%
%   yq = newtonDividedDiffInterp(x, y, xq)
%
%   Builds the same divided-difference table as the labsheet's Newton
%   polynomial algorithm,
%       d(k,1) = f(k)
%       d(k,j) = (d(k,j-1) - d(k-1,j-1)) / (x(k) - x(k-j+1)),  j = 2..k
%   then evaluates
%       P(x) = d(1,1) + d(2,2)(x-x1) + d(3,3)(x-x1)(x-x2) + ...
%   at every point in xq using nested multiplication (no explicit
%   product loop needed at evaluation time).
%
%   INPUTS
%     x  - [n x 1] or [1 x n] distinct node locations
%     y  - [n x 1] or [1 x n] f(x) at those nodes
%     xq - scalar or array of query points
%
%   OUTPUT
%     yq - interpolated value(s), same size as xq

    x = x(:);
    y = y(:);
    n = numel(x);
    if numel(y) ~= n
        error('newtonDividedDiffInterp:sizeMismatch', 'x and y must have the same length.');
    end

    % ---- Build the divided-difference table (only the diagonal is kept) ----
    D = zeros(n, n);
    D(:, 1) = y;
    for j = 2:n
        for k = j:n
            D(k, j) = (D(k, j-1) - D(k-1, j-1)) / (x(k) - x(k-j+1));
        end
    end
    coeffs = diag(D);   % [d(1,1); d(2,2); ...; d(n,n)]

    % ---- Nested evaluation of the Newton form at every query point ----
    yq = coeffs(n) * ones(size(xq));
    for k = n-1:-1:1
        yq = coeffs(k) + (xq - x(k)) .* yq;
    end
end
