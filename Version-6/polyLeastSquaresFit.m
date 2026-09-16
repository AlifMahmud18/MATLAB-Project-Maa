function coeffs = polyLeastSquaresFit(x, y, degree)
%POLYLEASTSQUARESFIT Least-squares polynomial regression via the normal
%   equations (EEE 212 Exp. 4), solved with gaussJordanSolve.m (Exp. 5)
%   instead of MATLAB's built-in polyfit/backslash.
%
%   coeffs = polyLeastSquaresFit(x, y, degree)
%
%   For degree = 1 this reduces to exactly the labsheet's linear
%   regression normal equations
%       n*a0 + (sum xi)*a1       = sum yi
%       (sum xi)*a0 + (sum xi^2)*a1 = sum xi*yi
%   and for general degree d it is the direct extension shown for the
%   quadratic case:
%       sum_{k=0}^{d} (sum xi^(j+k)) * a_k = sum xi^j * yi ,  j = 0..d
%
%   INPUTS
%     x, y   - paired data vectors
%     degree - polynomial degree d (1 = line, 2 = quadratic, ...)
%
%   OUTPUT
%     coeffs - [ (d+1) x 1 ] coefficients in ASCENDING power order,
%              i.e. y_m = coeffs(1) + coeffs(2)*x + coeffs(3)*x^2 + ...
%              (opposite convention to MATLAB's polyval, which expects
%              descending order - flip with coeffs(end:-1:1) if needed)

    x = x(:);
    y = y(:);
    n = numel(x);
    d = degree;

    % Power sums up to x^(2d), reused across the normal-equation matrix
    powerSums = zeros(1, 2*d + 1);
    for p = 0:2*d
        powerSums(p+1) = sum(x.^p);
    end

    M = zeros(d+1, d+1);
    rhs = zeros(d+1, 1);
    for j = 0:d
        for k = 0:d
            M(j+1, k+1) = powerSums(j+k+1);
        end
        rhs(j+1) = sum((x.^j) .* y);
    end

    coeffs = gaussJordanSolve(M, rhs);
end
