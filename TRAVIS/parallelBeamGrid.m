function [startRPZGrid, endRPZGrid] = parallelBeamGrid(startRPZ, endRPZ, n, m, spacing, varargin)
%PARALLELBEAMGRID Create an n-by-m grid of beams around an initial beam.
% Beams can be parallel or tilted toward the center line ("focus") by a given angle.
%
% Inputs (positional):
%   startRPZ  - 1x3 vector [R, phi, Z] for the start point (phi in radians)
%   endRPZ    - 1x3 vector [R, phi, Z] for the end point (phi in radians)
%   n, m      - integers, grid dimensions (rows n, columns m)
%   spacing   - scalar, Euclidean distance between adjacent beams (in Cartesian space)
%
% Name-Value parameters:
%   'RefAxis'        - 1x3 vector defining a reference axis to orient the grid (default [0 0 1])
%   'FocusAngleDeg'  - scalar angle in degrees to tilt each off-center beam toward the center line (default 0)
%                      For each beam offset by vector d in the plane perpendicular to the original beam,
%                      the direction becomes cos(a)*vhat - sin(a)*dhat, where a is the angle and dhat = d/|d|.
%                      The original beam (d = 0) remains unchanged.
%   'DoPlot'         - logical, if true plot the original and grid beams (default false)
%   'Axes'           - axes handle to plot into (default: create new figure)
%   'StartColor'     - RGB color for the original beam (default [0 0.4470 0.7410])
%   'GridColor'      - RGB color for the grid beams (default [0.8500 0.3250 0.0980])
%   'LineWidth'      - scalar line width for plots (default 1.5)
%
% Outputs:
%   startRPZGrid - n x m x 3 array of start points in [R, phi, Z]
%   endRPZGrid   - n x m x 3 array of end points in [R, phi, Z]
%
% Notes:
% - The grid lies in the plane perpendicular to the initial beam direction.
% - Distances are computed in Cartesian space; results are converted back to [R, phi, Z].
% - If n and m are odd, the center element equals the original beam location.
% - 'FocusAngleDeg' changes direction only (not length); beam length equals the original.

    % Parse name-value inputs
    p = inputParser;
    p.addParameter('RefAxis', [0 0 1], @(x) validateattributes(x, {'numeric'}, {'vector','numel',3}));
    p.addParameter('FocusAngleDeg', 0, @(x) validateattributes(x, {'numeric'}, {'scalar','nonnegative'}));
    p.addParameter('DoPlot', false, @(x) islogical(x) || isnumeric(x));
    p.addParameter('Axes', [], @(x) isempty(x) || ishghandle(x, 'axes'));
    p.addParameter('StartColor', [0 0.4470 0.7410], @(x) isnumeric(x) && numel(x)==3);
    p.addParameter('GridColor', [0.8500 0.3250 0.0980], @(x) isnumeric(x) && numel(x)==3);
    p.addParameter('LineWidth', 1.5, @(x) validateattributes(x, {'numeric'}, {'scalar','positive'}));
    p.parse(varargin{:});
    refAxis       = p.Results.RefAxis;
    focusAngleDeg = p.Results.FocusAngleDeg;
    doPlot        = logical(p.Results.DoPlot);
    ax            = p.Results.Axes;
    startColor    = p.Results.StartColor;
    gridColor     = p.Results.GridColor;
    lineWidth     = p.Results.LineWidth;

    % Convert start and end to Cartesian
    sXYZ = rpz2xyz(startRPZ);
    eXYZ = rpz2xyz(endRPZ);

    % Beam direction and length
    v = eXYZ - sXYZ;
    vnorm = norm(v);
    if vnorm == 0
        error('Start and end points coincide; beam direction is undefined.');
    end
    vhat = v / vnorm;

    % Orthonormal basis for the perpendicular plane
    r0 = refAxis(:).';
    if norm(r0) == 0
        r0 = [0 0 1];
    end
    r0 = r0 / norm(r0);
    if abs(dot(vhat, r0)) > 0.99
        r0 = [1 0 0];
    end
    u1 = cross(vhat, r0);
    u1n = norm(u1);
    if u1n < 1e-12
        r0 = [0 1 0];
        u1 = cross(vhat, r0);
        u1n = norm(u1);
        if u1n < 1e-12
            error('Failed to construct perpendicular basis.');
        end
    end
    u1 = u1 / u1n;
    u2 = cross(vhat, u1);
    u2 = u2 / norm(u2);

    % Precompute offsets
    rowOffsets = ((0:(n-1)) - (n-1)/2);  % size n
    colOffsets = ((0:(m-1)) - (m-1)/2);  % size m

    % Prepare outputs
    startRPZGrid = zeros(n, m, 3);
    endRPZGrid   = zeros(n, m, 3);

    % Store XYZ for plotting if needed
    sXYZGrid = zeros(n, m, 3);
    eXYZGrid = zeros(n, m, 3);

    % Angle in radians for focusing
    a = deg2rad(focusAngleDeg);
    if focusAngleDeg >= 90
        warning('FocusAngleDeg >= 90 degrees will reverse/side-aim beams; consider using a smaller angle.');
    end

    % Build grid
    for i = 1:n
        for j = 1:m
            d = spacing * (rowOffsets(i) * u1 + colOffsets(j) * u2);  % offset in perpendicular plane
            sShift = sXYZ + d;

            if focusAngleDeg > 0 && norm(d) > 0
                dhat = d / norm(d);
                w = cos(a) * vhat - sin(a) * dhat;  % tilt toward center line
                w = w / norm(w);
                eShift = sShift + vnorm * w;        % preserve original length
            else
                eShift = sShift + v;                % parallel case
            end

            % Store XYZ (for plotting)
            sXYZGrid(i,j,:) = sShift;
            eXYZGrid(i,j,:) = eShift;

            % Convert back to RPZ
            startRPZGrid(i,j,:) = xyz2rpz(sShift);
            endRPZGrid(i,j,:)   = xyz2rpz(eShift);
        end
    end

    % Optional plotting
    if doPlot
        if isempty(ax)
            fig = figure('Color','w');
            ax = axes('Parent',fig);
        end
        plotBeams(ax, sXYZ, eXYZ, sXYZGrid, eXYZGrid, startColor, gridColor, lineWidth, focusAngleDeg);
    end
end

function xyz = rpz2xyz(rpz)
    R = rpz(1);
    phi = rpz(2);
    Z = rpz(3);
    x = R * cos(phi);
    y = R * sin(phi);
    z = Z;
    xyz = [x, y, z];
end

function rpz = xyz2rpz(xyz)
    x = xyz(1); y = xyz(2); z = xyz(3);
    R = hypot(x, y);
    phi = atan2(y, x); % radians in [-pi, pi]
    Z = z;
    rpz = [R, phi, Z];
end

function plotBeams(ax, sXYZ, eXYZ, sXYZGrid, eXYZGrid, startColor, gridColor, lineWidth, focusAngleDeg)
    hold(ax, 'on');
    % Original beam
    plot3(ax, [sXYZ(1) eXYZ(1)], [sXYZ(2) eXYZ(2)], [sXYZ(3) eXYZ(3)], ...
          '-', 'Color', startColor, 'LineWidth', max(lineWidth, 2));

    % Grid beams
    [n, m, ~] = size(sXYZGrid);
    for i = 1:n
        for j = 1:m
            p1 = squeeze(sXYZGrid(i,j,:)).';
            p2 = squeeze(eXYZGrid(i,j,:)).';
            plot3(ax, [p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], ...
                  '-', 'Color', gridColor, 'LineWidth', lineWidth);
        end
    end

    % Mark the central start point (if exists in grid)
    ci = round((n+1)/2);
    cj = round((m+1)/2);
    p0 = squeeze(sXYZGrid(ci,cj,:)).';
    scatter3(ax, p0(1), p0(2), p0(3), 50, startColor, 'filled');

    axis(ax, 'equal');
    grid(ax, 'on');
    xlabel(ax, 'X'); ylabel(ax, 'Y'); zlabel(ax, 'Z');
    title(ax, sprintf('Beam Grid (FocusAngle = %.1f°)', focusAngleDeg));
    legend(ax, {'Original beam', 'Grid beams'}, 'Location', 'best');
    hold(ax, 'off');
end