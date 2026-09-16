function animateVibration(pos, bonds, modeShape, omega, amplitude, duration, fps)
%ANIMATEVIBRATION Play a real-time animation of one vibrational normal
%   mode of the crystal.
%
%   animateVibration(pos, bonds, modeShape, omega, amplitude, duration, fps)
%
%   PHYSICS BACKGROUND
%   In the harmonic approximation, once you know a mode's shape (the
%   displacement direction of each atom) and its angular frequency omega,
%   every atom's actual position over time is:
%
%       position(t) = equilibrium_position + amplitude * sin(omega*t) * modeShape
%
%   All atoms swing in and out of phase together (or exactly opposite,
%   depending on sign) at the SAME frequency omega - that's the defining
%   property of a normal mode.
%
%   INPUTS
%     pos        - [Natoms x 3] equilibrium atom positions
%     bonds      - [Nbonds x 5] bond list (from generateLattice)
%     modeShape  - [Natoms x 3] displacement pattern (from extractModeShape)
%     omega      - scalar angular frequency of this mode
%     amplitude  - (optional) how far atoms swing, in the same length
%                  units as pos. Default 0.25.
%     duration   - (optional) total animation length in seconds. Default 8.
%     fps        - (optional) frames per second. Default 30.
%
%   PERFORMANCE NOTE: rather than calling plot3/scatter3 fresh every
%   frame (slow - it destroys and recreates graphics objects), we create
%   each line/point object ONCE before the loop, then just update its
%   XData/YData/ZData each frame. This is the standard MATLAB technique
%   for smooth animation.

    if nargin < 5 || isempty(amplitude), amplitude = 0.25; end
    if nargin < 6 || isempty(duration),  duration  = 8;    end
    if nargin < 7 || isempty(fps),       fps       = 30;   end

    Natoms = size(pos, 1);
    Nbonds = size(bonds, 1);

    % --- Fix the axis limits ahead of time so the plot doesn't jitter/
    %     rescale every frame as atoms move ---
    maxAtomSwing = amplitude * max(sqrt(sum(modeShape.^2, 2)));
    margin = maxAtomSwing * 1.3 + 0.15;
    lo = min(pos, [], 1) - margin;
    hi = max(pos, [], 1) + margin;

    fig = figure('Name', 'Vibrating Crystal', 'Color', 'w');
    ax = axes(fig);
    hold(ax, 'on');
    axis(ax, 'equal');
    xlim(ax, [lo(1) hi(1)]);
    ylim(ax, [lo(2) hi(2)]);
    zlim(ax, [lo(3) hi(3)]);
    grid(ax, 'on');
    view(ax, 3);
    xlabel(ax, 'x'); ylabel(ax, 'y'); zlabel(ax, 'z');
    title(ax, sprintf('Normal mode animation (\\omega = %.4f)', omega));

    % --- Create graphics handles ONCE ---
    bondHandles = gobjects(Nbonds, 1);
    for k = 1:Nbonds
        i = bonds(k, 1);
        j = bonds(k, 2);
        bondHandles(k) = plot3(ax, ...
            [pos(i,1) pos(j,1)], [pos(i,2) pos(j,2)], [pos(i,3) pos(j,3)], ...
            'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
    end
    atomHandle = scatter3(ax, pos(:,1), pos(:,2), pos(:,3), 140, ...
        [0.85 0.33 0.10], 'filled', 'MarkerEdgeColor', 'k');

    % --- Animation loop ---
    nFrames = round(duration * fps);
    dt = 1 / fps;

    for f = 1:nFrames
        t = f * dt;
        displacedPos = pos + amplitude * sin(omega * t) * modeShape;

        set(atomHandle, ...
            'XData', displacedPos(:,1), ...
            'YData', displacedPos(:,2), ...
            'ZData', displacedPos(:,3));

        for k = 1:Nbonds
            i = bonds(k, 1);
            j = bonds(k, 2);
            set(bondHandles(k), ...
                'XData', [displacedPos(i,1) displacedPos(j,1)], ...
                'YData', [displacedPos(i,2) displacedPos(j,2)], ...
                'ZData', [displacedPos(i,3) displacedPos(j,3)]);
        end

        drawnow;
        pause(dt);

        % Allow early exit if the user closes the figure mid-animation
        if ~isvalid(fig)
            break;
        end
    end
end