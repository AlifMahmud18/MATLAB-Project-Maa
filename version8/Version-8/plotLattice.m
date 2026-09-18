function plotLattice(pos, bonds, ax)
%PLOTLATTICE Draw atoms (points) and bonds (lines) in 3D: plotLattice(pos, bonds, ax).
    if nargin < 3
        ax = gca;
    end
    hold(ax, 'on');
    for k = 1:size(bonds, 1)
        i = bonds(k, 1);
        j = bonds(k, 2);
        plot3(ax, [pos(i,1) pos(j,1)], [pos(i,2) pos(j,2)], [pos(i,3) pos(j,3)], ...
            'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
    end
    scatter3(ax, pos(:,1), pos(:,2), pos(:,3), 120, [0.85 0.33 0.10], ...
        'filled', 'MarkerEdgeColor', 'k');
    if isa(ax, 'matlab.ui.control.UIAxes')
        ax.DataAspectRatio = [1 1 1];
    else
        axis(ax, 'equal');
    end
    grid(ax, 'on');
    xlabel(ax, 'x'); ylabel(ax, 'y'); zlabel(ax, 'z');
    view(ax, 3);
    hold(ax, 'off');
end
