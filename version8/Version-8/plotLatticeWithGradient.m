function plotLatticeWithGradient(ax, coords, elementTypes, force_dist, bonds, cbLabel)
    %PLOTLATTICEWITHGRADIENT Draw the (deformed) lattice with atoms coloured by force_dist (jet colormap + one reused colorbar).
    %  plotLatticeWithGradient(ax, coords, elementTypes, force_dist, bonds)
    %  bonds is optional but recommended: bond lines between the deformed coords make the shear visible.

    hold(ax, 'on');

    if nargin >= 5 && ~isempty(bonds)
        for b = 1:size(bonds, 1)
            i = bonds(b, 1);
            j = bonds(b, 2);
            plot3(ax, [coords(i,1) coords(j,1)], [coords(i,2) coords(j,2)], [coords(i,3) coords(j,3)], ...
                'Color', [0.6 0.6 0.6], 'LineWidth', 1.0);
        end
    end

    scatter3(ax, coords(:,1), coords(:,2), coords(:,3), 150, force_dist, 'filled', 'MarkerEdgeColor', 'k');

    colormap(ax, jet);  clim(ax, [0 1]);
    if nargin < 6 || isempty(cbLabel), cbLabel = 'Relative Force Magnitude'; end
    cb = getappdata(ax, 'gradColorbar');
    if isempty(cb) || ~isvalid(cb), cb = colorbar(ax); setappdata(ax, 'gradColorbar', cb); end
    cb.Label.String = cbLabel;

    grid(ax, 'on');
    axis(ax, 'equal');
    view(ax, 3);
    hold(ax, 'off');
end
