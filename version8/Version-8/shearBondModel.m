function [f, kb, E] = shearBondModel(delta, k, model, alphaMorse)
%SHEARBONDMODEL Bond force f, stiffness kb = V'' and energy E for stretch delta = d - d0 (tension positive).
%  [f, kb, E] = shearBondModel(delta, k, model, alphaMorse)
%  'harmonic': V = k*delta^2/2.   'morse': V = D*(1 - exp(-a*delta))^2 with D = k/(2a^2).
%  Morse stiffness falls with stretch and reaches zero at delta = ln2/a (force maximum = bond rupture).

    switch lower(model)
        case 'harmonic'
            f  = k * delta;
            kb = k * ones(size(delta));
            E  = 0.5 * k * delta.^2;
        % ASSUMED bond law: Morse shape with D = k/(2a^2); calibrate k and alphaMorse to the material (bond energy, bond length, ideal strain).
        case 'morse'
            a  = alphaMorse;
            ex = exp(-a .* delta);
            f  = (k ./ a) .* ex .* (1 - ex);
            kb = k * (2 * ex.^2 - ex);
            E  = (k ./ (2 * a.^2)) .* (1 - ex).^2;
        otherwise
            error('shearBondModel:badModel', 'model must be ''harmonic'' or ''morse''.');
    end
end
