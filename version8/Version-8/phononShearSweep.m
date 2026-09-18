function res = phononShearSweep(pos, bonds, masses, opts)
%PHONONSHEARSWEEP Sweep xz shear (u_x = gamma*z, as in applyShearForce) and track phonons and bond breaking.
%  res = phononShearSweep(pos, bonds, masses, opts)   opts: k, gammaMax, nSteps, model, prestress, allowBreak, relax, alphaMorse.
%  Per strain: rotate bonds, relax internal atoms, remove Morse bonds past their force maximum,
%  build K -> D, project out the 3 translations, eigen-solve; lambda = omega^2 (< 0 = unstable).
%  res holds lambda, lowest eigenvectors, tau_xz, bond stretch/active flags, first-break and instability strains, omegaSoft (lowest 5% of modes, RMS).

    if nargin < 4, opts = struct(); end
    o = fillDefaults(opts);

    N  = size(pos, 1);
    nB = size(bonds, 1);
    gam = linspace(0, o.gammaMax, o.nSteps)';
    nS = numel(gam);

    d0 = bonds(:, 3);
    r0 = d0 .* bonds(:, 4:6);
    S  = [0 0 1; 0 0 0; 0 0 0];

    aMorse = o.alphaMorse ./ d0;
    if strcmpi(o.model, 'morse')
        deltaCrit = log(2) ./ aMorse;
    else
        deltaCrit = o.breakStretch * d0;
    end

    mv   = kron(masses(:), [1; 1; 1]);
    Minv = 1 ./ sqrt(mv);
    T = zeros(3*N, 3);
    for c = 1:3, T(c:3:end, c) = sqrt(masses(:)); end
    T = T ./ sqrt(sum(T.^2, 1));
    [Qf, ~] = qr(T);
    Q = Qf(:, 4:end);
    nModes = 3*N - 3;
    % Analysis choice, not a measured value: the soft band is the lowest 5% of modes.
    nSoft  = max(3, ceil(0.05 * nModes));

    Vol = N * (abs(det(o.latticeVectors)) / size(o.basisFrac, 1));

    res.gamma = gam;
    res.lambda = zeros(nModes, nS);
    res.modes = zeros(3*N, o.nLow, nS, 'single');
    res.tau = zeros(nS, 1);
    res.energy = zeros(nS, 1);
    res.stretch = zeros(nB, nS);
    res.damaged = false(nB, nS);
    res.active = true(nB, nS);
    res.rvec = zeros(nB, 3, nS);
    res.w = zeros(N, 3, nS, 'single');
    res.nFloppy = zeros(nS, 1);
    res.nUnstable = zeros(nS, 1);
    res.relaxResidual = zeros(nS, 1);
    res.omegaSoft = zeros(nS, 1);

    active  = true(nB, 1);
    damaged = false(nB, 1);
    w = zeros(N, 3);

    for s = 1:nS
        rAff = r0 + gam(s) * (r0 * S');
        resid = 0;
        for pass = 1:o.maxBreakPasses
            if o.relax
                [w, resid] = relaxNonaffine(w, rAff, bonds, d0, active, o, aMorse, N);
            end
            r = rAff + w(bonds(:, 2), :) - w(bonds(:, 1), :);
            d = sqrt(sum(r.^2, 2));
            delta = d - d0;
            over = delta > deltaCrit;
            damaged = damaged | over;
            newBreaks = over & active;
            if o.allowBreak && any(newBreaks)
                active(newBreaks) = false;
            else
                break
            end
        end

        [f, kb, E] = shearBondModel(delta, o.k, o.model, aMorse);
        f(~active) = 0;  kb(~active) = 0;  E(~active) = 0;

        res.tau(s)    = sum(f .* r(:, 1) .* r(:, 3) ./ d) / Vol;
        res.energy(s) = sum(E);
        res.stretch(:, s) = delta ./ d0;
        res.damaged(:, s) = damaged;
        res.active(:, s)  = active;
        res.rvec(:, :, s) = r;
        res.w(:, :, s)    = single(w);
        res.relaxResidual(s) = resid;

        if nModes == 0, continue; end
        K = assembleShearHessian(bonds, r, kb, f, N, o.prestress);
        D = K .* (Minv * Minv');
        Dp = Q' * D * Q;
        Dp = (Dp + Dp') / 2;
        if o.useJacobi
            [Vp, lam] = jacobiEigenSolver(Dp);
        else
            [Vp, L] = eig(Dp);
            lam = diag(L);
        end
        [lam, ord] = sort(lam);  Vp = Vp(:, ord);
        res.lambda(:, s) = lam;
        nl = min(o.nLow, nModes);
        res.modes(:, 1:nl, s) = single(Q * Vp(:, 1:nl));
        zt = o.zeroTolRel * max(abs(lam));
        res.nFloppy(s)   = sum(abs(lam) < zt);
        res.nUnstable(s) = sum(lam < -zt);
        ms = mean(lam(1:nSoft));
        res.omegaSoft(s) = sign(ms) * sqrt(abs(ms));
    end

    firstBreak = nan(nB, 1);
    for b = 1:nB
        idx = find(res.damaged(b, :), 1);
        if ~isempty(idx), firstBreak(b) = gam(idx); end
    end
    res.gammaFirstBreak = firstBreak;
    iu = find(res.nUnstable > 0, 1);
    if isempty(iu), res.gammaInstab = NaN; else, res.gammaInstab = gam(iu); end

    res.pos = pos;  res.bonds = bonds;  res.masses = masses;
    res.opts = o;   res.deltaCrit = deltaCrit;  res.aMorse = aMorse;
    res.N = N;      res.Minv = Minv;  res.nSoft = nSoft;
end

function o = fillDefaults(opts)
    % ASSUMED - calibrate to the real material: alphaMorse (critical stretch = ln2/alpha = 17.3%) and breakStretch (harmonic model);
    % fit them to the material's ideal tensile strain / bond energy. They set every break and Morse-instability strain.
    o = struct('k', 1, 'gammaMax', 0.5, 'nSteps', 40, 'model', 'morse', ...
        'prestress', true, 'allowBreak', true, 'alphaMorse', 4, ...
        'breakStretch', 0.15, 'relax', true, 'nLow', 6, 'useJacobi', false, ...
        'zeroTolRel', 1e-8, 'maxBreakPasses', 8, 'relaxIter', 60, 'relaxTol', 1e-9, ...
        'latticeVectors', eye(3), 'basisFrac', [0 0 0]);
    fn = fieldnames(opts);
    for i = 1:numel(fn), o.(fn{i}) = opts.(fn{i}); end
end

function [w, resid] = relaxNonaffine(w, rAff, bonds, d0, active, o, aMorse, N)
    i = bonds(:, 1);  j = bonds(:, 2);
    dmax = 0.1 * min(d0);
    for it = 1:o.relaxIter
        [F, kb, f, r, E0] = forcesAt(w);
        resid = max(abs(F(:)));
        if resid < o.relaxTol, return; end
        K = assembleShearHessian(bonds, r, kb, f, N, true);
        [V, L] = eig(K);
        lam = diag(L);
        keep = abs(lam) > 1e-8 * max(abs(lam));
        Fv = reshape(F', [], 1);
        dw = V(:, keep) * ((V(:, keep)' * Fv) ./ abs(lam(keep)));
        step = reshape(dw, 3, N)';
        mx = max(sqrt(sum(step.^2, 2)));
        if mx > dmax, step = step * (dmax / mx); end
        slope = -sum(F(:) .* step(:));
        t = 1;
        while t > 1e-8
            E1 = energyAt(w + t * step);
            if E1 <= E0 + 1e-4 * t * slope, break; end
            t = t / 2;
        end
        w = w + t * step;
    end
    [F, ~] = forcesAt(w);
    resid = max(abs(F(:)));

    function [F, kb, f, r, E] = forcesAt(wv)
        r = rAff + wv(j, :) - wv(i, :);
        d = sqrt(sum(r.^2, 2));
        [f, kb, Eb] = shearBondModel(d - d0, o.k, o.model, aMorse);
        f(~active) = 0;  kb(~active) = 0;  Eb(~active) = 0;
        E = sum(Eb);
        fv = (f ./ d) .* r;
        F = zeros(N, 3);
        for c = 1:3
            F(:, c) = accumarray(i, fv(:, c), [N 1]) - accumarray(j, fv(:, c), [N 1]);
        end
    end
    function E = energyAt(wv)
        d = sqrt(sum((rAff + wv(j, :) - wv(i, :)).^2, 2));
        [~, ~, Eb] = shearBondModel(d - d0, o.k, o.model, aMorse);
        Eb(~active) = 0;
        E = sum(Eb);
    end
end
