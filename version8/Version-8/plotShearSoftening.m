function figs = plotShearSoftening(S, outDir)
%PLOTSHEARSOFTENING Six phonon-softening figures P1-P6 from phononShearSuite results.
%  P1 soft-band curve, P2 DOS shift, P3 spectral map, P4 softest-mode atoms, P5 stress vs phonon, P6 transport.
%  Data: eigenvalues/eigenvectors from phononShearSweep, extractModeShape, centralDiffDerivative, phononTransportMetrics.
%  Curves stop at the first strain with omega^2 < 0; later points describe a different structure.
%  plotShearSoftening(S, outDir) also saves PNGs to outDir if given.

    if nargin < 2, outDir = ''; end
    if ~isempty(outDir) && ~exist(outDir, 'dir'), mkdir(outDir); end
    C = S.cases;  M = C{S.main};  g = M.gamma;  nS = numel(g);
    col = [0.00 0.45 0.74; 0.47 0.67 0.19; 0.85 0.33 0.10; 0.49 0.18 0.56; 0.64 0.08 0.18];
    fs = 12;  lw = 2.2;

    gInst = M.gammaInstab;  gBreak = min(M.gammaFirstBreak);
    if isnan(gInst), nV = nS; else, nV = max(find(g < gInst, 1, 'last'), 2); end
    sShow = nV;  gShow = g(sShow);
    valid = 1:nV;
    omega = sqrt(max(M.lambda, 0));
    wSoft0 = M.omegaSoft(1);

    figs.P1 = newFig([100 100 1150 470]);
    subplot(1, 2, 1); hold on;
    plotValid(g, M.omegaSoft / wSoft0, nV, col(3, :), lw, 'soft band (lowest 5% of modes, RMS)');
    wm = mean(omega, 1)' / mean(omega(:, 1));
    plot(g(valid), wm(valid), '-', 'Color', col(1, :), 'LineWidth', 1.6, 'DisplayName', 'mean of all modes');
    events(gBreak, gInst); xlabel('Shear strain \gamma'); ylabel('\omega / \omega(0)');
    title(sprintf('Softening under the Morse model (\\alpha = %g, critical stretch %.1f%%)', M.opts.alphaMorse, 100 * log(2) / M.opts.alphaMorse)); legend('Location', 'southwest'); styleAx(fs);
    subplot(1, 2, 2); hold on;
    for c = 1:numel(C)
        Cc = C{c}; ni = nValidOf(Cc);
        plotValid(Cc.gamma, Cc.omegaSoft / Cc.omegaSoft(1), ni, col(c, :), lw, S.labels{c});
    end
    events(gBreak, gInst); xlabel('Shear strain \gamma'); ylabel('soft-band \omega / \omega(0)');
    title('Same result across three model levels'); legend('Location', 'southwest', 'FontSize', fs - 2); styleAx(fs);
    finish(figs.P1, outDir, 'P1_softening_curve');

    figs.P2 = newFig([100 100 900 480]); hold on;
    wmax = max(omega(:, valid), [], 'all');
    wg = linspace(0, 1.08 * wmax, 700)';  sig = 0.015 * wmax;
    d0 = gaussDOS(omega(:, 1), wg, sig);  dS = gaussDOS(omega(:, sShow), wg, sig);
    area(wg, d0, 'FaceColor', col(1, :), 'FaceAlpha', 0.25, 'EdgeColor', col(1, :), 'LineWidth', lw, 'DisplayName', '\gamma = 0');
    area(wg, dS, 'FaceColor', col(3, :), 'FaceAlpha', 0.25, 'EdgeColor', col(3, :), 'LineWidth', lw, ...
        'DisplayName', sprintf('\\gamma = %.3f (last stable)', gShow));
    m0 = mean(omega(:, 1));  mS = mean(omega(:, sShow));
    xline(m0, '--', 'Color', col(1, :), 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xline(mS, '--', 'Color', col(3, :), 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlabel('\omega'); ylabel('g(\omega)');
    title(sprintf('Phonon density of states  (mean \\omega %+.1f%%, lowest 5%% of modes %+.1f%%)', ...
        100 * (mS - m0) / m0, 100 * (M.omegaSoft(sShow) - wSoft0) / wSoft0));
    legend('Location', 'northeast'); styleAx(fs);
    finish(figs.P2, outDir, 'P2_dos_shift');

    figs.P3 = newFig([100 100 950 500]);
    map = zeros(numel(wg), nV);
    for s = 1:nV, map(:, s) = gaussDOS(omega(:, s), wg, sig); end
    imagesc(g(valid), wg, map); axis xy; hold on;
    colormap(gca, flipud(gray)); cb = colorbar; cb.Label.String = 'g(\omega)';
    plot(g(valid), omega(1, valid), 'r-', 'LineWidth', lw);
    xlabel('Shear strain \gamma'); ylabel('\omega (red = lowest mode)');
    title('Spectral map: every phonon band vs shear (stable range only)'); styleAx(fs);
    finish(figs.P3, outDir, 'P3_spectral_map');

    figs.P4 = newFig([100 100 1150 500]);
    for q = 1:2
        s = [1, sShow];  s = s(q);
        subplot(1, 2, q); hold on;
        pos = M.pos;  N = M.N;
        cur = pos + g(s) * (pos * [0 0 1; 0 0 0; 0 0 0]') + double(M.w(:, :, s));
        u = extractModeShape(double(M.modes(:, 1, s)), M.masses, N);
        mag = sqrt(sum(u.^2, 2));
        [~, ord] = sort(mag);
        scatter(cur(ord, 1), cur(ord, 3), 40 + 330 * mag(ord), mag(ord), 'filled', 'MarkerEdgeColor', [0.2 0.2 0.2]);
        quiver(cur(:, 1), cur(:, 3), u(:, 1), u(:, 3), 0.5, 'Color', 'k', 'LineWidth', 0.9, 'MaxHeadSize', 0.5);
        colormap(gca, turbo); clim([0 1]); axis equal; box on;
        xlabel('x'); ylabel('z');
        title(sprintf('\\gamma = %.3f,  \\omega_1 = %.4g', g(s), omega(1, s)));
        cb = colorbar; cb.Label.String = '|u_i| / max|u|'; styleAx(fs);
    end
    sgtitle('Softest vibrational mode: which atoms move (all atoms projected on xz, arrows = in-plane displacement)', 'FontSize', fs, 'FontWeight', 'bold');
    finish(figs.P4, outDir, 'P4_softest_mode_atoms');

    [Gt, G0] = tangentModulus(g, M.tau);
    figs.P5 = newFig([100 100 1000 500]);
    yyaxis left; hold on;
    plot(g(valid), M.tau(valid), '-', 'LineWidth', lw, 'Color', col(1, :));
    ylabel('Macroscopic shear stress \tau_{xz}'); ax = gca; ax.YColor = col(1, :);
    yyaxis right; hold on;
    plot(g(valid), M.omegaSoft(valid) / wSoft0, '-', 'LineWidth', lw, 'Color', col(3, :));
    plot(g(valid), sqrt(max(Gt(valid) / G0, 0)), '--', 'LineWidth', lw, 'Color', col(4, :));
    ylabel('\omega_{soft}/\omega_{soft}(0)   and   v_s/v_{s0} = \surd(G_t/G_0)'); ax.YColor = [0.3 0.3 0.3];
    legend({'\tau_{xz}', 'phonon soft band', 'shear-wave speed \surd(G_t/G_0)'}, 'Location', 'southwest');
    events(gBreak, gInst); xlabel('Shear strain \gamma');
    title('Mechanical response vs phonon softening'); styleAx(fs);
    finish(figs.P5, outDir, 'P5_stress_vs_phonon');

    % ASSUMED transport inputs (15.5 THz, 77/300/450 K): replace with the material's real top phonon frequency and the device temperatures.
    tm = phononTransportMetrics(M.lambda(:, valid), 15.5, [77 300 450]);
    gv = g(valid);
    figs.P6 = newFig([100 100 1400 450]);
    tcol = [0.00 0.45 0.74; 0.47 0.67 0.19; 0.85 0.33 0.10];
    subplot(1, 3, 1); hold on;
    for t = 1:3, plot(gv, tm.u2Ratio(:, t), '-', 'LineWidth', lw, 'Color', tcol(t, :), 'DisplayName', sprintf('%g K', tm.Tlist(t))); end
    set(gca, 'YScale', 'log'); xlabel('Shear strain \gamma'); ylabel('\langle u^2\rangle / \langle u^2\rangle_0');
    title('Atomic mean-square displacement'); legend('Location', 'northwest'); styleAx(fs);
    subplot(1, 3, 2); hold on;
    for t = 1:3, plot(gv, tm.muDW(:, t), '-', 'LineWidth', lw, 'Color', tcol(t, :), 'DisplayName', sprintf('\\mu_{DW}, %g K', tm.Tlist(t))); end
    plot(gv, max(Gt(valid) / G0, 0), 'k--', 'LineWidth', lw, 'DisplayName', '\mu_{ac} \propto G_t/G_0');
    xlabel('Shear strain \gamma'); ylabel('\mu / \mu_0 (phonon-limited)'); ylim([0 max(1.1, max(tm.muDW(:)))]);
    title('Relative phonon-limited mobility'); legend('Location', 'southwest'); styleAx(fs);
    subplot(1, 3, 3); hold on;
    plot(gv, tm.vSat, '-', 'LineWidth', lw, 'Color', col(4, :), 'DisplayName', 'v_{sat}/v_{sat,0}');
    plot(gv, tm.thetaDRatio, '-', 'LineWidth', lw, 'Color', col(2, :), 'DisplayName', '\Theta_D^{eff}/\Theta_{D,0}');
    xlabel('Shear strain \gamma'); ylabel('ratio to unstrained'); title('Saturation velocity & Debye temperature proxies');
    legend('Location', 'best'); styleAx(fs);
    sgtitle('Consequences for carrier transport (ratios to unstrained crystal, stable strain range)', 'FontSize', fs, 'FontWeight', 'bold');
    finish(figs.P6, outDir, 'P6_transport');
end

function nV = nValidOf(M)
    if isnan(M.gammaInstab), nV = numel(M.gamma); else, nV = max(find(M.gamma < M.gammaInstab, 1, 'last'), 2); end
end

function plotValid(g, y, nV, c, lw, name)
    plot(g(1:nV), y(1:nV), '-', 'Color', c, 'LineWidth', lw, 'DisplayName', name);
    if nV < numel(g)
        plot(g(nV:min(nV + 1, end)), y(nV:min(nV + 1, end)), ':', 'Color', c, 'LineWidth', lw, 'HandleVisibility', 'off');
        xline(g(nV + 1), '-', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
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

function events(gBreak, gInst)
    yl = ylim;
    if ~isnan(gBreak)
        xline(gBreak, ':', sprintf('first bond breaks %.3f', gBreak), 'Color', [0.85 0.33 0.10], 'LineWidth', 1.6, ...
            'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    end
    if ~isnan(gInst)
        xline(gInst, '--', sprintf('\\omega^2<0 at %.3f', gInst), 'Color', [0.64 0.08 0.18], 'LineWidth', 1.6, ...
            'LabelVerticalAlignment', 'middle', 'HandleVisibility', 'off');
    end
    ylim(yl);
end

function y = gaussDOS(v, grid, sig)
    y = zeros(size(grid));
    for q = 1:numel(v), y = y + exp(-0.5 * ((grid - v(q)) / sig).^2); end
    y = y / (sig * sqrt(2 * pi));
end

function f = newFig(pos), f = figure('Color', 'w', 'Position', pos); end
function styleAx(fs), set(gca, 'FontSize', fs, 'LineWidth', 1.1, 'Box', 'on'); grid on; end

function finish(f, outDir, name)
    drawnow;
    set(findall(f, 'Type', 'axes'), 'Toolbar', []);
    if ~isempty(outDir), exportgraphics(f, fullfile(outDir, [name '.png']), 'Resolution', 170); end
end
