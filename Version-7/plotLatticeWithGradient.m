function plotLatticeWithGradient(ax, coords, elementTypes, force_dist)
    % PLOTLATTICEWITHGRADIENT Renders atoms colored by local force distribution gradient
    
    % Ensure axes holds the plot
    hold(ax, 'on');
    
    % Scatter plot mapping force_dist to the colormap
    scatter3(ax, coords(:,1), coords(:,2), coords(:,3), 150, force_dist, 'filled', 'MarkerEdgeColor', 'k');
    
    % Apply standard JET colormap and add colorbar
    colormap(ax, jet);
    cb = colorbar(ax);
    cb.Label.String = 'Relative Force Magnitude';
    
    % Formatting
    grid(ax, 'on');
    axis(ax, 'equal');
    view(ax, 3);
    hold(ax, 'off');
end