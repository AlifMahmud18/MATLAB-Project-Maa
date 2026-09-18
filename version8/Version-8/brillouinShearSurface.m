function fig = brillouinShearSurface(latticeVectors, basisFrac, basisMasses, cutoff, k, gammaMax, alphaMorse)
%BRILLOUINSHEARSURFACE 3D phonon dispersion surface w(kx,ky) of a sheared crystal, with a shear-strain slider.
%   brillouinShearSurface(latticeVectors, basisFrac, basisMasses, cutoff, k, gammaMax, alphaMorse)  (no arguments = simple cubic demo)
%   Bonds: every neighbour pair of the unit cell within cutoff, images included; D(k) = M^-1/2 K(k) M^-1/2, K(k) has phases exp(i k.r).
%   Strain steps come from phononShearSweep (bond rotation, relaxation, Morse rupture); slider value = engineering shear gamma.
%   Surface = lowest branch(es) on the kz = 0 slice, masked to the first Brillouin zone of the sheared lattice; w < 0 = imaginary.
%   Failure: an instability at small |k| = long sliding planes, at the zone edge = atom-by-atom alternating shear.

    if nargin < 1 || isempty(latticeVectors)
        [latticeVectors, basisFrac, basisMasses] = getBuiltinCrystalDef('sc');  cutoff = 1.5;
    end
    if nargin < 5 || isempty(k), k = 1; end
    if nargin < 6 || isempty(gammaMax), gammaMax = 0.5; end
    if nargin < 7 || isempty(alphaMorse), alphaMorse = 4; end

    LV = latticeVectors;  nb = size(basisFrac, 1);  masses = basisMasses(:);
    bonds = cellNeighborBonds(LV, basisFrac, cutoff);
    nB = size(bonds, 1);  d0 = bonds(:, 3);  bi = bonds(:, 1);  bj = bonds(:, 2);
    posCell = basisFrac * LV;
    Minv = kron(1 ./ sqrt(masses), [1; 1; 1]);
    kmax = 0.55 * max(vecnorm(2 * pi * inv(LV)', 2, 2));
    dim = 3 * nb;
    if dim <= 30, nk = 41; elseif dim <= 72, nk = 31; else, nk = 23; end
    kv = linspace(-kmax, kmax, nk);  [KX, KY] = meshgrid(kv, kv);
    Kall = [KX(:), KY(:), zeros(numel(KX), 1)];
    [IA, IB, JA, JB] = deal(zeros(nB, 9));
    for a = 1:3, for b = 1:3
        c = 3 * (a - 1) + b;  IA(:, c) = 3*(bi-1) + a;  IB(:, c) = 3*(bi-1) + b;
        JA(:, c) = 3*(bj-1) + a;  JB(:, c) = 3*(bj-1) + b;
    end, end
    rowsT = [IA(:); JA(:); IA(:); JA(:)];  colsT = [IB(:); JB(:); JB(:); IB(:)];

    nSteps = 101;  gam = linspace(0, gammaMax, nSteps)';
    modelNames = {'Harmonic, affine (slider model)', 'Harmonic + pre-stress + relaxation', 'Morse + pre-stress + relaxation + rupture'};
    modelOpts = {{'harmonic', false, false, false}, {'harmonic', true, false, true}, {'morse', true, true, true}};
    sweeps = cell(1, 3);  refLam = cell(1, 3);  refMask = cell(1, 3);  cacheLam = cell(3, nSteps);
    mIdx = 3;  busy = false;  floppy = false;

    fig = figure('Name', 'Brillouin-zone dispersion surface under shear', 'NumberTitle', 'off', 'Position', [100 60 1280 740], 'Color', 'w');
    ax = axes('Parent', fig, 'Units', 'pixels', 'Position', [75 60 800 640]);
    pnl = uipanel(fig, 'Units', 'pixels', 'Position', [960 15 305 710], 'Title', 'Shear control', 'BackgroundColor', 'w');
    mk = @(style, str, pos, varargin) uicontrol(pnl, 'Style', style, 'String', str, 'Position', pos, 'BackgroundColor', 'w', varargin{:});
    mk('text', 'Engineering shear strain gamma', [15 670 275 20], 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
    sld = mk('slider', '', [15 648 275 18], 'Min', 0, 'Max', gammaMax, 'Value', 0, 'SliderStep', [0.01 0.1], 'Callback', @(s, ~) onSlide(s.Value));
    addlistener(sld, 'ContinuousValueChange', @(s, ~) onSlide(s.Value));
    gLabel = mk('text', 'gamma = 0.0000', [15 622 275 20], 'HorizontalAlignment', 'left');
    mk('text', 'Bond model', [15 594 275 18], 'HorizontalAlignment', 'left');
    mk('popupmenu', modelNames, [15 568 275 24], 'Value', 3, 'Callback', @(s, ~) onModel(s.Value));
    cbAll = mk('checkbox', 'Show the 3 lowest branches', [15 540 275 22], 'Value', 0, 'Callback', @(~, ~) onSlide(sld.Value));
    mk('pushbutton', 'Jump to first zero-touch', [15 500 275 32], 'FontWeight', 'bold', 'ForegroundColor', 'w', ...
        'BackgroundColor', [0.75 0.15 0.15], 'Callback', @(~, ~) jumpToZeroTouch());
    mk('pushbutton', 'Reset view', [15 465 275 28], 'Callback', @(~, ~) view(ax, -35, 28));
    txt = mk('text', '', [15 15 275 440], 'HorizontalAlignment', 'left', 'FontName', 'Consolas', 'FontSize', 9);
    rotate3d(fig, 'on');
    ensureModel(mIdx);
    lam0 = refLam{mIdx};  mask0 = refMask{mIdx};  floppy = isFloppy(lam0, mask0);
    wTop = 1.15 * max(sqrt(max(lam0(:, 3), 0)));  w1max = max(sqrt(max(lam0(:, 1), 0)));
    zlo = -0.6 * w1max;  zhi = wTop;
    hSurf = gobjects(1, 3);
    Z0 = nan(nk);
    for q = 1:3
        hSurf(q) = surf(ax, KX, KY, Z0, Z0, 'EdgeColor', [0.25 0.25 0.25], 'EdgeAlpha', 0.2, 'FaceAlpha', 1 - 0.55 * (q > 1));
        hold(ax, 'on');
    end
    surf(ax, [-kmax kmax; -kmax kmax], [-kmax -kmax; kmax kmax], zeros(2), 'FaceColor', [0.4 0.4 0.4], 'FaceAlpha', 0.18, 'EdgeColor', 'none');
    hZero = scatter3(ax, NaN, NaN, NaN, 9, 'r', 'filled');
    hMin = scatter3(ax, NaN, NaN, NaN, 160, 'd', 'filled', 'MarkerFaceColor', [1 0.9 0], 'MarkerEdgeColor', 'k', 'LineWidth', 1.4);
    nneg = round(256 * (-zlo) / (zhi - zlo));
    colormap(ax, [[linspace(0.35, 1, nneg)', zeros(nneg, 1), zeros(nneg, 1)]; parula(256 - nneg)]);
    clim(ax, [zlo zhi]);  zlim(ax, [zlo zhi]);  xlim(ax, [-kmax kmax]);  ylim(ax, [-kmax kmax]);
    xlabel(ax, 'k_x');  ylabel(ax, 'k_y');  zlabel(ax, '\omega  (imaginary shown negative)');
    view(ax, -35, 28);  grid(ax, 'on');  cb = colorbar(ax);  cb.Label.String = '\omega';
    onSlide(0);

    function onSlide(g)
        if busy, return; end
        busy = true;  cleaner = onCleanup(@() setBusy(false));
        [~, s] = min(abs(gam - g));  update(s);
    end
    function setBusy(v), busy = v; end

    function onModel(m)
        mIdx = m;  ensureModel(m);  lam0 = refLam{m};  mask0 = refMask{m};  floppy = isFloppy(lam0, mask0);  onSlide(sld.Value);
    end

    function ensureModel(m)
        if ~isempty(sweeps{m}), return; end
        txt.String = {'Computing the strained cell...'};  drawnow;
        o = struct('k', k, 'gammaMax', gammaMax, 'nSteps', nSteps, 'model', modelOpts{m}{1}, 'prestress', modelOpts{m}{2}, ...
            'allowBreak', modelOpts{m}{3}, 'relax', modelOpts{m}{4}, 'alphaMorse', alphaMorse, 'latticeVectors', LV, 'basisFrac', basisFrac);
        sweeps{m} = phononShearSweep(posCell, bonds, masses, o);
        refMask{m} = bzMetric(0, Kall) <= 1 + 1e-9;
        lam = nan(size(Kall, 1), 3);  lam(refMask{m}, :) = dispersion(m, 1, Kall(refMask{m}, :));
        refLam{m} = lam;
    end
    function [lam, mf] = dispersion(m, s, Kin, variant)
        if nargin < 4, variant = 'full'; end
        res = sweeps{m};  r = res.rvec(:, :, s);  d = sqrt(sum(r.^2, 2));  u = r ./ d;
        delta = res.stretch(:, s) .* d0;
        [f, kb] = shearBondModel(delta, k, modelOpts{m}{1}, res.aMorse);
        act = res.active(:, s);  if strcmp(variant, 'noRupture'), act(:) = true; end
        f(~act) = 0;  kb(~act) = 0;
        g = f ./ d;  if ~modelOpts{m}{2}, g(:) = 0; end
        switch variant
            case 'noCompPre', g(g < 0) = 0;
            case 'noStrPre', g(g > 0) = 0;
            case 'noPre', g(:) = 0;
            case 'noSoften', kb(delta > 0 & act) = k;
        end
        B = zeros(nB, 9);
        for a = 1:3, for b = 1:3
            B(:, 3*(a-1)+b) = (kb - g) .* u(:, a) .* u(:, b) + g * (a == b);
        end, end
        nKp = size(Kin, 1);  lam = zeros(nKp, 3);
        linK = rowsT + (colsT - 1) * dim;  nEnt = numel(linK);  chunk = 200;
        for i0 = 1:chunk:nKp
            idxK = i0:min(i0 + chunk - 1, nKp);  nC = numel(idxK);
            phK = reshape(exp(1i * (r * Kin(idxK, :)')), nB, 1, nC);
            Bp = reshape(B .* phK, nB * 9, nC);  Bc = reshape(B .* conj(phK), nB * 9, nC);
            valsK = [repmat(B(:), 1, nC); repmat(B(:), 1, nC); -Bp; -Bc];
            Hs = sparse(repmat(linK, nC, 1), repelem((1:nC)', nEnt, 1), valsK(:), dim * dim, nC);
            Hp = reshape(full(Hs), dim, dim, nC) .* (Minv * Minv');
            Hp = (Hp + pagectranspose(Hp)) / 2;
            evK = sort(real(reshape(pageeig(Hp, 'vector'), dim, nC)), 1);
            lam(idxK, :) = evK(1:3, :)';
        end
        mf = [];
    end

    function mf = bzMetric(gs, Kin)
        LVs = LV;  LVs(:, 1) = LVs(:, 1) + gs * LV(:, 3);
        Bs = 2 * pi * inv(LVs)';
        [n1, n2, n3] = ndgrid(-2:2);  n = [n1(:) n2(:) n3(:)];  n(all(n == 0, 2), :) = [];
        G = n * Bs;  mf = max((Kin * G') ./ (sum(G.^2, 2)' / 2), [], 2);
    end

    function update(s)
        gs = gam(s);  res = sweeps{mIdx};  gLabel.String = sprintf('gamma = %.4f', gs);
        mf = bzMetric(gs, Kall);  inside = mf <= 1 + 1e-9;
        lam = cacheLam{mIdx, s};
        if isempty(lam)
            lam = nan(numel(mf), 3);  lam(inside, :) = dispersion(mIdx, s, Kall(inside, :));  cacheLam{mIdx, s} = lam;
        end
        showAll = cbAll.Value;
        for q = 1:3
            W = sign(lam(:, q)) .* sqrt(abs(lam(:, q)));  W = reshape(W, nk, nk);
            set(hSurf(q), 'ZData', W, 'CData', W, 'Visible', matlabOnOff(q == 1 || showAll));
        end
        tol = 1e-9 * max(lam0(:, 1));
        kn = vecnorm(Kall, 2, 2);  notG = kn > 1e-9 & inside;
        neg = notG & lam(:, 1) < -tol;
        w1 = sign(lam(:, 1)) .* sqrt(abs(lam(:, 1)));
        set(hZero, 'XData', Kall(neg, 1), 'YData', Kall(neg, 2), 'ZData', w1(neg));
        cmp = notG & mask0 & lam0(:, 1) > 1e-8 * max(lam0(:, 1));
        ratio = nan(size(mf));  ratio(cmp) = lam(cmp, 1) ./ lam0(cmp, 1);
        [rMin, iR] = min(ratio);
        if isempty(iR) || isnan(rMin)
            set(hMin, 'XData', NaN, 'YData', NaN, 'ZData', NaN);
        else
            set(hMin, 'XData', Kall(iR, 1), 'YData', Kall(iR, 2), 'ZData', w1(iR));
        end
        [lMin, iL] = min(lam(notG, 1));  idxNG = find(notG);  iL = idxNG(iL);
        nBroken = sum(~res.active(:, s));
        unstable = any(neg);
        lines = {sprintf('gamma = %.4f   (%s)', gs, shortName(mIdx)), sprintf('Morse alpha = %g (critical stretch %.1f%%)', alphaMorse, 100 * log(2) / alphaMorse), sprintf('broken bonds: %d of %d', nBroken, nB), ''};
        if ~isnan(rMin)
            lines = [lines, {sprintf('softest wave (lowest branch):'), sprintf('  k = (%+.3f, %+.3f)', Kall(iR, 1), Kall(iR, 2)), ...
                sprintf('  %.0f%% of the way to the zone edge', 100 * mf(iR)), sprintf('  w^2 / w0^2 = %.3f', rMin), ...
                sprintf('lowest w^2 anywhere = %+.4g', lMin), sprintf('  at k = (%+.3f, %+.3f)', Kall(iL, 1), Kall(iL, 2)), ''}];
        end
        if isempty(iR) || isnan(rMin), frR = NaN; else, frR = mf(iR); end
        if unstable
            fr = frR;
            lines = [lines, {'STATUS: UNSTABLE - surface has', 'crossed the zero floor.', sprintf('%d k-points with w^2 < 0 (red dots)', nnz(neg)), ...
                failureText(fr)}];
            if ~isempty(iL), lines = [lines, {''}, explainInstability(s, Kall(iL, :))]; end
        else
            lines = [lines, {'STATUS: stable (surface above zero)', 'yellow diamond = softest wave.', 'It becomes the failure mode when it', 'touches the grey zero plane.', ...
                failureText(frR)}];
        end
        if floppy, lines = [lines, {'', 'NOTE: this spring network is nearly floppy', 'along some k directions, so small', 'strains can already push w^2 below 0.'}]; end
        txt.String = lines;
        title(ax, sprintf('\\gamma = %.4f    min \\omega^2 (k\\neq0) = %+.4g', gs, lMin));
        drawnow limitrate;
    end

    function jumpToZeroTouch()
        txt.String = {'Scanning strains for the first zero-touch...'};  drawnow;
        kc = linspace(-kmax, kmax, 17);  [cx, cy] = meshgrid(kc, kc);  Kc = [cx(:), cy(:), zeros(numel(cx), 1)];
        found = 0;
        for s = 1:3:nSteps
            if unstableAt(s, Kc), found = s;  break; end
        end
        if found == 0 && unstableAt(nSteps, Kc), found = nSteps; end
        if found == 0
            txt.String = [{'No zero-touch up to gamma_max.'}, txt.String(:)'];  return
        end
        s = max(found - 2, 1);
        while s < found && ~unstableAt(s, Kc), s = s + 1; end
        sld.Value = gam(s);  busy = false;  onSlide(gam(s));
    end

    function tf = unstableAt(s, Kc)
        mf = bzMetric(gam(s), Kc);  ins = mf <= 1 + 1e-9 & vecnorm(Kc, 2, 2) > 1e-9;
        lm = dispersion(mIdx, s, Kc(ins, :));
        tf = any(lm(:, 1) < -1e-9 * max(lam0(:, 1)));
    end

    function lines = explainInstability(s, kpt)
        variants = {'noRupture', 'noCompPre', 'noStrPre', 'noSoften', 'noPre'};
        tt = dispersion(mIdx, s, kpt, 'full');  lFull = tt(1);
        lv = zeros(1, 5);
        for iv = 1:5, tt = dispersion(mIdx, s, kpt, variants{iv});  lv(iv) = tt(1); end
        dl = lv - lFull;  small = 0.05 * abs(lFull);
        res = sweeps{mIdx};  st = res.stretch(:, s);  act = res.active(:, s);
        comp = st < -0.02 & act;  strt = st > 0.02 & act;
        if strcmp(modelOpts{mIdx}{1}, 'morse'), crit = 100 * log(2) / alphaMorse; else, crit = 15; end
        detail = {sprintf('loss of %d ruptured bonds', nnz(~act)), ...
            sprintf('compressive pre-stress of %d squeezed bonds (mean %.1f%%): buckling-like', nnz(comp), 100 * mean(st(comp))), ...
            sprintf('pre-stress of %d stretched bonds', nnz(strt)), ...
            sprintf('softening of %d stretched bonds (mean +%.1f%%, critical %.1f%%)', nnz(strt), 100 * mean(st(strt)), crit), ...
            'the combined pre-stress of all loaded bonds'};
        dmax = max(dl(1:4));
        if dmax > small
            idx = find(dl(1:4) >= 0.5 * dmax & dl(1:4) > small);
        elseif dl(5) > small
            idx = 5;
        else
            lines = {'REASON: no tested ingredient (rupture,', 'pre-stress, Morse softening) matters:', 'the bond network stiffness itself has', 'run out (geometry / connectivity).'};
            return
        end
        [~, o] = sort(dl(idx), 'descend');  idx = idx(o);
        lines = {'REASON (largest effect first):'};
        for iv = idx
            if lv(iv) >= 0, tag = ' -> restores stability'; else, tag = ' -> less unstable'; end
            lines{end+1} = ['  - ' detail{iv} tag]; %#ok<AGROW>
        end
        if numel(idx) > 1 && all(lv(idx) >= 0), lines{end+1} = '  (removing any one of these is enough)'; end
    end
    function n = shortName(m), n = strtok(modelNames{m}, ','); end
    function tf = isFloppy(l0, m0)
        kk = vecnorm(Kall, 2, 2);  ok = m0 & kk > 1e-9;  c = l0(ok, 1) ./ kk(ok).^2;
        % Heuristic: below 2% of the median acoustic stiffness the network is called nearly floppy.
        tf = min(c) < 0.02 * median(c);
    end
end

function s = failureText(fr)
    % Heuristic cut-offs (35% / 75% of the way to the zone edge) for the failure signature; not calibrated to a material.
    if isnan(fr)
        s = '';
    elseif fr < 0.35
        s = 'Failure signature: long wavelength -> macroscopic sliding planes.';
    elseif fr < 0.75
        s = 'Failure signature: intermediate wavelength -> mesoscopic domains / twins.';
    else
        s = 'Failure signature: zone-edge wave -> atom-by-atom alternating shear.';
    end
end

function s = matlabOnOff(tf)
    if tf, s = 'on'; else, s = 'off'; end
end
