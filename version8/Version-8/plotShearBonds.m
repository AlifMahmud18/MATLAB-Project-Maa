function figs = plotShearBonds(S, outDir, parent)
%PLOTSHEARBONDS Bond-breaking plot groups Q1-Q2 from phononShearSuite results.
%  Q1: all bonds projected on xz at three strains, coloured by stretch, broken bonds black (where).
%  Q2: broken % per bond length vs strain (when) and stretch vs n_x*n_z at the last stable strain (why).
%  To first order stretch = gamma*n_x*n_z, so bonds on the 45-degree tension axis break first.
%  Only strains before the first omega^2 < 0 are drawn. Works for any crystal.
%  plotShearBonds(S) opens Q1-Q2 as two figures.
%  plotShearBonds(S, outDir) also saves PNGs to outDir.
%  plotShearBonds(S, outDir, parent) opens no window at all: it fills `parent` (a uitab,
%  uipanel or uifigure) with a tab group holding Q1-Q2, replacing whatever was in there.
%  Returns a struct of the containers (figure or tab) keyed Q1, Q2.

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

    M = S.cases{S.main};  g = M.gamma;  b = M.bonds;  nB = size(b, 1);
    fs = 12;  lw = 2;
    gInst = M.gammaInstab;  gFirst = min(M.gammaFirstBreak);
    nS = numel(g);
    if isnan(gInst), nV = nS; else, nV = max(find(g < gInst, 1, 'last'), 2); end

    shellVal = unique(round(b(:, 3), 2));
    shellId = arrayfun(@(d) find(abs(shellVal - d) < 0.006, 1), b(:, 3));
    shellCol = turbo(max(numel(shellVal), 2) + 2);  shellCol = shellCol(2:end-1, :);

    noBreak = isnan(gFirst) || gFirst > g(nV);
    if noBreak
        snap = unique(max(2, round(nV * [0.5 0.8 1])));
        while numel(snap) < 3, snap(end + 1) = min(nV, snap(end) + 1); end
    else
        s1 = max(2, find(g >= 0.5 * gFirst, 1));
        s2 = min(nV, find(g >= gFirst, 1));
        snap = unique([s1, s2, nV]);
        while numel(snap) < 3, snap(end + 1) = min(nV, snap(end) + 1); end
        snap = snap(1:3);
    end

    ttl = sprintf('Where shear breaks bonds (xz projection; Morse \\alpha = %g, critical stretch %.1f%%)', M.opts.alphaMorse, 100 * log(2) / M.opts.alphaMorse);
    if noBreak
        ttl = sprintf(['%s' newline 'No bond reaches its critical stretch (max %.1f%% vs %.1f%%) before the lattice goes unstable (omega^2<0 at gamma = %.3f): elastic failure, not rupture'], ...
            ttl, 100 * max(max(M.stretch(:, 1:nV))), 100 * log(2) / M.opts.alphaMorse, gInst);
    end

    cv = shearPlotCanvas(tg, 'Q1 Bond map', 1, 3, [80 80 1450 480], ttl);
    figs.Q1 = cv.container;
    smax = max(0.05, max(max(M.stretch(:, 1:nV))));  clim0 = [-smax, smax];
    cm = turbo(256);  S0 = [0 0 1; 0 0 0; 0 0 0];
    for q = 1:3
        s = snap(q);  ax = cv.ax(q);  hold(ax, 'on');
        cur = M.pos + g(s) * (M.pos * S0') + double(M.w(:, :, s));
        r = M.rvec(:, :, s);
        [~, order] = sort(abs(M.stretch(:, s)));
        for k = order'
            P = cur(b(k, 1), [1 3]);  Qp = P + r(k, [1 3]);
            if M.active(k, s)
                ci = round(1 + 255 * min(max((M.stretch(k, s) - clim0(1)) / (clim0(2) - clim0(1)), 0), 1));
                plot(ax, [P(1) Qp(1)], [P(2) Qp(2)], '-', 'Color', [cm(ci, :) 0.75], 'LineWidth', 1.2);
            else
                plot(ax, [P(1) Qp(1)], [P(2) Qp(2)], '--', 'Color', 'k', 'LineWidth', 2);
                mid = 0.5 * (P + Qp);  plot(ax, mid(1), mid(2), 'kx', 'MarkerSize', 8, 'LineWidth', 1.8);
            end
        end
        plot(ax, cur(:, 1), cur(:, 3), 'o', 'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerEdgeColor', 'w', 'MarkerSize', 6);
        colormap(ax, turbo); clim(ax, clim0); ax.DataAspectRatio = [1 1 1]; box(ax, 'on');
        xlabel(ax, 'x'); ylabel(ax, 'z');
        nBr = sum(~M.active(:, s));
        title(ax, sprintf('\\gamma = %.3f   (%d broken bonds, %.1f%%)', g(s), nBr, 100 * nBr / nB));
        set(ax, 'FontSize', fs, 'LineWidth', 1.1);
        if q == 3, cb = colorbar(ax); cb.Label.String = 'bond stretch (d-d_0)/d_0    [black x = broken]'; end
    end
    cv.finish(outDir, 'Q1_bond_map');

    cv = shearPlotCanvas(tg, 'Q2 Break statistics', 1, 2, [80 80 1150 470]);
    figs.Q2 = cv.container;
    ax = cv.ax(1); hold(ax, 'on');
    broken = ~M.active;
    plot(ax, g(1:nV), 100 * mean(broken(:, 1:nV), 1), 'k-', 'LineWidth', 2.6, 'DisplayName', 'all bonds');
    for sh = 1:numel(shellVal)
        sel = shellId == sh;
        plot(ax, g(1:nV), 100 * mean(broken(sel, 1:nV), 1), '-', 'LineWidth', lw, 'Color', shellCol(sh, :), ...
            'DisplayName', sprintf('d_0 = %.2f (%d bonds)', shellVal(sh), nnz(sel)));
    end
    xlabel(ax, 'Shear strain \gamma'); ylabel(ax, 'Broken bonds (%)'); title(ax, 'When: bond-breaking progression by bond length');
    lg = legend(ax, 'Location', 'northwest'); lg.FontSize = fs - 3;
    set(ax, 'FontSize', fs, 'LineWidth', 1.1); grid(ax, 'on'); box(ax, 'on');

    ax = cv.ax(2); hold(ax, 'on');
    s = nV;
    n = b(:, 4:6);  geo = n(:, 1) .* n(:, 3);
    for sh = 1:numel(shellVal)
        sel = shellId == sh & M.active(:, s);
        scatter(ax, geo(sel), 100 * M.stretch(sel, s), 36, shellCol(sh, :), 'filled', 'MarkerFaceAlpha', 0.7);
    end
    brk = ~M.active(:, s);
    scatter(ax, geo(brk), 100 * M.stretch(brk, s), 90, 'k', 'x', 'LineWidth', 1.8);
    xx = linspace(min(geo), max(geo), 20);
    plot(ax, xx, 100 * g(s) * xx, 'k--', 'LineWidth', 1.4);
    if strcmp(M.opts.model, 'morse')
        yline(ax, 100 * log(2) / M.opts.alphaMorse, ':', 'critical stretch', 'Color', [0.64 0.08 0.18], 'LineWidth', 1.6);
    end
    xlabel(ax, 'Geometric factor  n_x n_z  of the bond direction'); ylabel(ax, sprintf('Bond stretch at \\gamma = %.3f (%%)', g(s)));
    title(ax, 'Why: stretch \approx \gamma n_x n_z  (dashed, 1st order); x = broken');
    set(ax, 'FontSize', fs, 'LineWidth', 1.1); grid(ax, 'on'); box(ax, 'on');
    cv.finish(outDir, 'Q2_break_statistics');
end
