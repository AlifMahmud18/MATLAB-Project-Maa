function [I, Icum] = compositeTrapz(y, h)
%COMPOSITETRAPZ Composite trapezoidal rule on a uniform grid (EEE 212 Exp. 7).
%  [I, Icum] = compositeTrapz(y, h)   y: samples, h: spacing; I = total integral, Icum = running integral.

    y = y(:);
    n = numel(y);
    Icum = zeros(n, 1);
    for k = 2:n
        Icum(k) = Icum(k-1) + h * (y(k) + y(k-1)) / 2;
    end
    I = Icum(end);
end
