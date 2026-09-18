function idx = atomHover(ax, labels, posFcn)
%ATOMHOVER Show a label next to the atom under the mouse pointer (no click needed).
%  atomHover(ax, labels, posFcn): labels = {Natoms x 1} strings, posFcn() = current [Natoms x 3] positions.
%  The mouse ray (ax.CurrentPoint) is compared with every atom; the nearest one to the viewer
%  within 4.5% of the atom span gets the label. One WindowButtonMotionFcn serves all registered axes.
%  Call again after any cla() because the label is a child of the axes.

    idx = [];
    if ischar(ax) && strcmp(ax, 'pick')
        idx = pickAtom(labels, posFcn);
        return
    end
    fig = ancestor(ax, 'figure');
    lbl = text(ax, NaN, NaN, NaN, '', 'Visible', 'off', 'Clipping', 'off', ...
        'BackgroundColor', [1 1 0.85], 'EdgeColor', [0.2 0.2 0.2], 'Margin', 5, ...
        'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'none', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
        'HitTest', 'off', 'PickableParts', 'none');

    reg = getappdata(fig, 'atomHoverRegistry');
    if isempty(reg), reg = struct('ax', {}, 'lbl', {}, 'labels', {}, 'posFcn', {}); end
    keep = arrayfun(@(r) isvalid(r.ax) && r.ax ~= ax, reg);
    reg = reg(keep);
    reg(end + 1) = struct('ax', ax, 'lbl', lbl, 'labels', {labels}, 'posFcn', posFcn);
    setappdata(fig, 'atomHoverRegistry', reg);

    fig.WindowButtonMotionFcn = @(~, ~) onMove(fig);
end

function onMove(fig)
    reg = getappdata(fig, 'atomHoverRegistry');
    for r = reg
        if ~isvalid(r.ax) || ~isvalid(r.lbl), continue; end
        try
            P = r.posFcn();
            [i, span] = pickAtom(P, r.ax.CurrentPoint);
            if i == 0
                r.lbl.Visible = 'off';
            else
                r.lbl.String = r.labels{i};
                r.lbl.Position = P(i, :) + [0 0 0.03 * span];
                r.lbl.Visible = 'on';
            end
        catch
            r.lbl.Visible = 'off';
        end
    end
end

function [i, span] = pickAtom(P, cp)
    p1 = cp(1, :);  p2 = cp(2, :);
    dirv = (p2 - p1) / norm(p2 - p1);
    rel = P - p1;
    t = rel * dirv';
    perp = sqrt(max(sum(rel.^2, 2) - t.^2, 0));
    span = max(max(P) - min(P));
    if span <= 0, span = 1; end
    hit = find(perp < 0.045 * span);
    if isempty(hit)
        i = 0;
    else
        [~, k] = min(t(hit));
        i = hit(k);
    end
end
