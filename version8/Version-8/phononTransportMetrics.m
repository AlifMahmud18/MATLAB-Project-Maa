function tm = phononTransportMetrics(lambda, fRefTHz, Tlist)
%PHONONTRANSPORTMETRICS Turn a strained phonon spectrum into carrier-transport proxies (ratios to gamma = 0).
%  tm = phononTransportMetrics(lambda, fRefTHz, Tlist)   lambda: [nModes x nS] omega^2, column 1 = unstrained.
%  <u^2> ~ mean (1/omega)*coth(hbar*omega/2kT); mobility_DW = <u^2>_0/<u^2>; mobility_ac = (omega_low ratio)^2.
%  v_sat ~ sqrt(omega_op) (top 10% of modes); Debye ratio from the log-mean frequency; Bose occupation of the softest mode.
%  The top unstrained mode is mapped to fRefTHz (default 15.5); modes below 5% of the lowest are floored.

    % ASSUMED defaults: 15.5 THz top phonon (Si-like) and T = 77/300/450 K; use the real material's top phonon frequency and the device temperature.
    if nargin < 2 || isempty(fRefTHz), fRefTHz = 15.5; end
    if nargin < 3 || isempty(Tlist),   Tlist   = [77 300 450]; end

    hbar = 1.054571817e-34;  kB = 1.380649e-23;
    [nModes, nS] = size(lambda);
    omega = sqrt(max(lambda, 0));
    omega0max = max(omega(:, 1));
    % Chosen values: 5% floor on the lowest unstrained frequency; the top 10% of modes are treated as the optical-like end.
    omegaFloor = 0.05 * omega(1, 1);
    omegaEff = max(omega, omegaFloor);
    unstable = any(lambda < 1e-7, 1)';

    OmegaRef = 2 * pi * fRefTHz * 1e12;
    nT = numel(Tlist);
    u2 = zeros(nS, nT);  nBoseLow = zeros(nS, nT);
    for t = 1:nT
        x = hbar * OmegaRef * (omegaEff / omega0max) / (kB * Tlist(t));
        term = (1 ./ omegaEff) .* coth(x / 2);
        u2(:, t) = mean(term, 1)';
        nBoseLow(:, t) = 1 ./ (exp(x(1, :)') - 1);
    end
    tm.Tlist = Tlist;
    tm.u2Ratio       = u2 ./ u2(1, :);
    tm.rmsRatio      = sqrt(tm.u2Ratio);
    tm.muDW          = 1 ./ tm.u2Ratio;
    tm.muAc          = (omegaEff(1, :)' / omegaEff(1, 1)).^2;
    tm.nBoseLow      = nBoseLow;
    tm.nBoseLowRatio = nBoseLow ./ nBoseLow(1, :);

    nTop = max(1, ceil(0.10 * nModes));
    omegaOp = mean(omega(end-nTop+1:end, :), 1)';
    tm.vSat = sqrt(omegaOp / omegaOp(1));
    tm.omegaOpRatio = omegaOp / omegaOp(1);

    omegaLog = exp(mean(log(omegaEff), 1))';
    tm.thetaDRatio = omegaLog / omegaLog(1);

    tm.omegaLowRatio = omega(1, :)' / omega(1, 1);
    tm.unstable = unstable;
    tm.omega = omega;
    tm.fRefTHz = fRefTHz;
    tm.omega0max = omega0max;
    tm.hbarOmegaRef_meV = hbar * OmegaRef / 1.602176634e-19 * 1e3;
end
