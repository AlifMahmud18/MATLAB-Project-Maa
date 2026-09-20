function figs = plotShearSoftening(S, outDir, parent)
%PLOTSHEARSOFTENING Six phonon-softening plot groups P1-P6 from phononShearSuite results.
%  P1 soft-band curve, P2 DOS shift, P3 spectral map, P4 softest-mode atoms, P5 stress vs phonon, P6 transport.
%  Data: eigenvalues/eigenvectors from phononShearSweep, extractModeShape, centralDiffDerivative, phononTransportMetrics.
%  Curves stop at the first strain with omega^2 < 0; later points describe a different structure.
%  plotShearSoftening(S) opens P1-P6 as six figures.
%  plotShearSoftening(S, outDir) also saves PNGs to outDir.
%  plotShearSoftening(S, outDir, parent) opens no window at all: it fills `parent` (a uitab,
%  uipanel or uifigure) with a tab group holding P1-P6, replacing whatever was in there.
%  Returns a struct of the containers (figure or tab) keyed P1..P6.

    if nargin < 2, outDir = ''; end
    if nargin < 3, parent = []; end
    if ~isempty(outDir) && ~exist(outDir, 'dir'), mkdir(outDir); end

    tg = [];
    if ~isempty(parent)
        delete(allchild(parent));
        host = uigridlayout(parent, [1 1]);
        host.Padding = [0 0 0 0];
        tg = uitabgroup(host);
        tg.Layout.Row = 1;
        tg.Layout.Column = 1;
    end

    C = S.cases;  M = C{S.main};  g = M.gamma;  nS = numel(g);
    col = [0.00 0.45 0.74; 0.47 0.67 0.19; 0.85 0.33 0.10; 0.49 0.18 0.56; 0.64 0.08 0.18];
    fs = 12;  lw = 2.2;

    gInst = M.gammaInstab;  gBreak = min(M.gammaFirstBreak);
    if isnan(gInst), nV = nS; else, nV = max(find(g < gInst, 1, 'last'), 2); end
    sShow = nV;  gShow = g(sShow);
    valid = 1:nV;
    omega = sqrt(max(M.lambda, 0));
    wSoft0 = M.omegaSoft(1);

    cv = shearPlotCanvas(tg, 'P1 Softening curve', 1, 2, [100 100 1150 470]);
    figs.P1 = cv.container;
    ax = cv.ax(1);  hold(ax, 'on');
    plotValid(ax, g, M.omegaSoft / wSoft0, nV, col(3, :), lw, 'soft band (lowest 5% of modes, RMS)');
    wm = mean(omega, 1)' / mean(omega(:, 1));
    plot(ax, g(valid), wm(valid), '-', 'Color', col(1, :), 'LineWidth', 1.6, 'DisplayName', 'mean of all modes');
    markEvents(ax, gBreak, gInst); xlabel(ax, 'Shear strain \gamma'); ylabel(ax, '\omega / \omega(0)');
    title(ax, sprintf('Softening under the Morse model (\\alpha = %g, critical stretch %.1f%%)', M.opts.alphaMorse, 100 * log(2) / M.opts.alphaMorse));
    legend(ax, 'Location', 'southwest'); styleAx(ax, fs);
    ax = cv.ax(2);  hold(ax, 'on');
    for c = 1:numel(C)
        Cc = C{c}; ni = nValidOf(Cc);
        plotValid(ax, Cc.gamma, Cc.omegaSoft / Cc.omegaSoft(1), ni, col(c, :), lw, S.labels{c});
    end
    markEvents(ax, gBreak, gInst); xlabel(ax, 'Shear strain \gamma'); ylabel(ax, 'soft-band \omega / \omega(0)');
    title(ax, 'Same result across three model levels'); legend(ax, 'Location', 'southwest', 'FontSize', fs - 2); styleAx(ax, fs);
    cv.finish(outDir, 'P1_softening_curve');

    cv = shearPlotCanvas(tg, 'P2 DOS shift', 1, 1, [100 100 900 480]);
    figs.P2 = cv.container;
    ax = cv.ax(1);  hold(ax, 'on');
    wmax = max(omega(:, valid), [], 'all');
    wg = linspace(0, 1.08 * wmax, 700)';  sig = 0.015 * wmax;
    d0 = gaussDOS(omega(:, 1), wg, sig);  dS = gaussDOS(omega(:, sShow), wg, sig);
    area(ax, wg, d0, 'FaceColor', col(1, :), 'FaceAlpha', 0.25, 'EdgeColor', col(1, :), 'LineWidth', lw, 'DisplayName', '\gamma = 0');
    area(ax, wg, dS, 'FaceColor', col(3, :), 'FaceAlpha', 0.25, 'EdgeColor', col(3, :), 'LineWidth', lw, ...
        'DisplayName', sprintf('\\gamma = %.3f (last stable)', gShow));
    m0 = mean(omega(:, 1));  mS = mean(omega(:, sShow));
    xline(ax, m0, '--', 'Color', col(1, :), 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xline(ax, mS, '--', 'Color', col(3, :), 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlabel(ax, '\omega'); ylabel(ax, 'g(\omega)');
    title(ax, sprintf('Phonon density of states  (mean \\omega %+.1f%%, lowest 5%% of modes %+.1f%%)', ...
        100 * (mS - m0) / m0, 100 * (M.omegaSoft(sShow) - wSoft0) / wSoft0));
    legend(ax, 'Location', 'northeast'); styleAx(ax, fs);
    cv.finish(outDir, 'P2_dos_shift');

    cv = shearPlotCanvas(tg, 'P3 Spectral map', 1, 1, [100 100 950 500]);
    figs.P3 = cv.container;
    ax = cv.ax(1);
    map = zeros(numel(wg), nV);
    for s = 1:nV, map(:, s) = gaussDOS(omega(:, s), wg, sig); end
    imagesc(ax, g(valid), wg, map);  ax.YDir = 'normal';  hold(ax, 'on');
    colormap(ax, flipud(gray)); cb = colorbar(ax); cb.Label.String = 'g(\omega)';
    plot(ax, g(valid), omega(1, valid), 'r-', 'LineWidth', lw);
    xlabel(ax, 'Shear strain \gamma'); ylabel(ax, '\omega (red = lowest mode)');
    title(ax, 'Spectral map: every phonon band vs shear (stable range only)'); styleAx(ax, fs);
    cv.finish(outDir, 'P3_spectral_map');

    cv = shearPlotCanvas(tg, 'P4 Softest mode', 1, 2, [100 100 1150 500], ...
        'Softest vibrational mode: which atoms move (all atoms projected on xz, arrows = in-plane displacement)');
    figs.P4 = cv.container;
    for q = 1:2
        s = [1, sShow];  s = s(q);
        ax = cv.ax(q);  hold(ax, 'on');
        pos = M.pos;  N = M.N;
        cur = pos + g(s) * (pos * [0 0 1; 0 0 0; 0 0 0]') + double(M.w(:, :, s));
        u = extractModeShape(double(M.modes(:, 1, s)), M.masses, N);
        mag = sqrt(sum(u.^2, 2));
        [~, ord] = sort(mag);
        scatter(ax, cur(ord, 1), cur(ord, 3), 40 + 330 * mag(ord), mag(ord), 'filled', 'MarkerEdgeColor', [0.2 0.2 0.2]);
        quiver(ax, cur(:, 1), cur(:, 3), u(:, 1), u(:, 3), 0.5, 'Color', 'k', 'LineWidth', 0.9, 'MaxHeadSize', 0.5);
        colormap(ax, turbo); clim(ax, [0 1]); ax.DataAspectRatio = [1 1 1]; box(ax, 'on');
        xlabel(ax, 'x'); ylabel(ax, 'z');
        title(ax, sprintf('\\gamma = %.3f,  \\omega_1 = %.4g', g(s), omega(1, s)));
        cb = colorbar(ax); cb.Label.String = '|u_i| / max|u|'; styleAx(ax, fs);
    end
    cv.finish(outDir, 'P4_softest_mode_atoms');

    [Gt, G0] = tangentModulus(g, M.tau);
    cv = shearPlotCanvas(tg, 'P5 Stress vs phonon', 1, 1, [100 100 1000 500]);
    figs.P5 = cv.container;
    ax = cv.ax(1);
    yyaxis(ax, 'left'); hold(ax, 'on');
    plot(ax, g(valid), M.tau(valid), '-', 'LineWidth', lw, 'Color', col(1, :));
    ylabel(ax, 'Macroscopic shear stress \tau_{xz}'); ax.YColor = col(1, :);
    yyaxis(ax, 'right'); hold(ax, 'on');
    plot(ax, g(valid), M.omegaSoft(valid) / wSoft0, '-', 'LineWidth', lw, 'Color', col(3, :));
    plot(ax, g(valid), sqrt(max(Gt(valid) / G0, 0)), '--', 'LineWidth', lw, 'Color', col(4, :));
    ylabel(ax, '\omega_{soft}/\omega_{soft}(0)   and   v_s/v_{s0} = \surd(G_t/G_0)'); ax.YColor = [0.3 0.3 0.3];
    legend(ax, {'\tau_{xz}', 'phonon soft band', 'shear-wave speed \surd(G_t/G_0)'}, 'Location', 'southwest');
    markEvents(ax, gBreak, gInst); xlabel(ax, 'Shear strain \gamma');
    title(ax, 'Mechanical response vs phonon softening'); styleAx(ax, fs);
    cv.finish(outDir, 'P5_stress_vs_phonon');

    % ASSUMED transport inputs (15.5 THz, 77/300/450 K): replace with the material's real top phonon frequency and the device temperatures.
    tm = phononTransportMetrics(M.lambda(:, valid), 15.5, [77 300 450]);
    gv = g(valid);
    cv = shearPlotCanvas(tg, 'P6 Transport', 1, 3, [100 100 1400 450], ...
        'Consequences for carrier transport (ratios to unstrained crystal, stable strain range)');
    figs.P6 = cv.container;
    tcol = [0.00 0.45 0.74; 0.47 0.67 0.19; 0.85 0.33 0.10];
    ax = cv.ax(1); hold(ax, 'on');
    for t = 1:3, plot(ax, gv, tm.u2Ratio(:, t), '-', 'LineWidth', lw, 'Color', tcol(t, :), 'DisplayName', sprintf('%g K', tm.Tlist(t))); end
    set(ax, 'YScale', 'log'); xlabel(ax, 'Shear strain \gamma'); ylabel(ax, '\langle u^2\rangle / \langle u^2\rangle_0');
    title(ax, 'Atomic mean-square displacement'); legend(ax, 'Location', 'northwest'); styleAx(ax, fs);
    ax = cv.ax(2); hold(ax, 'on');
    for t = 1:3, plot(ax, gv, tm.muDW(:, t), '-', 'LineWidth', lw, 'Color', tcol(t, :), 'DisplayName', sprintf('\\mu_{DW}, %g K', tm.Tlist(t))); end
    plot(ax, gv, max(Gt(valid) / G0, 0), 'k--', 'LineWidth', lw, 'DisplayName', '\mu_{ac} \propto G_t/G_0');
    xlabel(ax, 'Shear strain \gamma'); ylabel(ax, '\mu / \mu_0 (phonon-limited)'); ylim(ax, [0 max(1.1, max(tm.muDW(:)))]);
    title(ax, 'Relative phonon-limited mobility'); legend(ax, 'Location', 'southwest'); styleAx(ax, fs);
    ax = cv.ax(3); hold(ax, 'on');
    plot(ax, gv, tm.vSat, '-', 'LineWidth', lw, 'Color', col(4, :), 'DisplayName', 'v_{sat}/v_{sat,0}');
    plot(ax, gv, tm.thetaDRatio, '-', 'LineWidth', lw, 'Color', col(2, :), 'DisplayName', '\Theta_D^{eff}/\Theta_{D,0}');
    xlabel(ax, 'Shear strain \gamma'); ylabel(ax, 'ratio to unstrained'); title(ax, 'Saturation velocity & Debye temperature proxies');
    legend(ax, 'Location', 'best'); styleAx(ax, fs);
    cv.finish(outDir, 'P6_transport');
end

function nV = nValidOf(M)
    if isnan(M.gammaInstab), nV = numel(M.gamma); else, nV = max(find(M.gamma < M.gammaInstab, 1, 'last'), 2); end
end

function plotValid(ax, g, y, nV, c, lw, name)
    plot(ax, g(1:nV), y(1:nV), '-', 'Color', c, 'LineWidth', lw, 'DisplayName', name);
    if nV < numel(g)
        plot(ax, g(nV:min(nV + 1, end)), y(nV:min(nV + 1, end)), ':', 'Color', c, 'LineWidth', lw, 'HandleVisibility', 'off');
        xline(ax, g(nV + 1), '-', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
    end
end

function [Gt, G0] = tangentModulus(g, tau)
    f = @(x) interp1(g, tau, min(max(x, g(1)), g(end)), 'pchip');
    h = 0.5 * (g(2) - g(1));
    Gt = centralDiffDerivative(f, g, h);
    Gt = Gt(:);
    Gt(1) = Gt(2);  Gt(end) = Gt(end-1);
    G0 = Gt(1);
end

function markEvents(ax, gBreak, gInst)
    yl = ylim(ax);
    if ~isnan(gBreak)
        xline(ax, gBreak, ':', sprintf('first bond breaks %.3f', gBreak), 'Color', [0.85 0.33 0.10], 'LineWidth', 1.6, ...
            'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    end
    if ~isnan(gInst)
        xline(ax, gInst, '--', sprintf('\\omega^2<0 at %.3f', gInst), 'Color', [0.64 0.08 0.18], 'LineWidth', 1.6, ...
            'LabelVerticalAlignment', 'middle', 'HandleVisibility', 'off');
    end
    ylim(ax, yl);
end

function y = gaussDOS(v, grid, sig)
    y = zeros(size(grid));
    for q = 1:numel(v), y = y + exp(-0.5 * ((grid - v(q)) / sig).^2); end
    y = y / (sig * sqrt(2 * pi));
end

function styleAx(ax, fs)
    set(ax, 'FontSize', fs, 'LineWidth', 1.1, 'Box', 'on');
    grid(ax, 'on');
end
