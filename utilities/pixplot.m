function hpcolor = pixplot(varargin)
% PIXPLOT Create a pseudocolor plot with pixel-centered ticks and no clipping.
%
% Usage:
%   hp = pixplot(V)
%   hp = pixplot(X, Y, V)
%   hp = pixplot(ax, V)
%   hp = pixplot(ax, X, Y, V)
%   hp = pixplot(..., 'Name', Value, ...)
%
% Inputs:
%   - ax: axes handle to plot into (optional; if omitted, uses current axes)
%   - V : numeric matrix [nx x ny] of values
%   - X : optional centers for x-dimension. Either:
%         * vector length nx (row or column), or
%         * matrix size [nx x ny] with columns identical (meshgrid-style)
%   - Y : optional centers for y-dimension. Either:
%         * vector length ny (row or column), or
%         * matrix size [nx x ny] with rows identical (meshgrid-style)
%
% Name-Value options:
%   'NumTicks'     : maximum number of ticks per axis (default 7)
%   'TickPrecision': number of decimals in tick labels (default 2)
%   'ShowColorbar' : true/false (default true)
%   'AxisImage'    : true/false to set axis image (default true)
%   'Colormap'     : colormap name or matrix (default [])
%   'CLim'         : [min max] color limits (default [])
%   'Mask'         : logical matrix same size as V; true = visible, false = transparent
%   'EdgeColor'    : color for patch edges (default 'none')
%   'XLabel'       : string label for x-axis (default '')
%   'YLabel'       : string label for y-axis (default '')
%   'Title'        : string title (default '')
%
% Output:
%   - hpcolor: handle to the pcolor patch object
%
% Notes:
%   - This function constructs edge coordinates from pixel centers to avoid
%     the clipping behavior of pcolor (which drops the last row/column).
%   - It assumes a Cartesian, axis-aligned grid (rectangular cells).
%   - If X/Y are matrices, they must be meshgrid-like (X varies only by row,
%     Y varies only by column). For general distorted grids, provide proper
%     vertex coordinates and use pcolor/surf directly.

    % ----------------------
    % Parse axes handle
    % ----------------------
    args = varargin;
    ax = [];
    if ~isnumeric(args{1})
        ax = args{1};
        args(1) = [];
    else
        ax = gca;
    end

    % ----------------------
    % Separate positional and name-value parts
    % ----------------------
    % Identify where name-value pairs start
    nvStart = find(cellfun(@(c) ischar(c) || isstring(c), args), 1, 'first');
    if isempty(nvStart)
        posArgs = args;
        nvArgs  = {};
    else
        posArgs = args(1:nvStart-1);
        nvArgs  = args(nvStart:end);
    end

    % ----------------------
    % Decode positional args: V or X,Y,V
    % ----------------------
    if numel(posArgs) == 1
        V = posArgs{1};
        validateattributes(V, {'numeric'}, {'2d'}, mfilename, 'V');
        [nx, ny] = size(V);
        xcenters = (1:nx).';
        ycenters = (1:ny);
    elseif numel(posArgs) == 3
        X = posArgs{1};
        Y = posArgs{2};
        V = posArgs{3};
        validateattributes(V, {'numeric'}, {'2d'}, mfilename, 'V');
        [nx, ny] = size(V);

        [xcenters, ycenters] = parseCentersXY(X, Y, nx, ny);
    else
        error('pixplot:InvalidInputs', 'Expected either V or X,Y,V (optionally with an axes handle).');
    end

    % ----------------------
    % Parse name-value options
    % ----------------------
    p = inputParser;
    p.addParameter('NumTicks', 7, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    p.addParameter('TickPrecision', 2, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    p.addParameter('ShowColorbar', true, @(b) islogical(b) || isnumeric(b));
    p.addParameter('AxisImage', true, @(b) islogical(b) || isnumeric(b));
    p.addParameter('Colormap', [], @(c) ischar(c) || isstring(c) || isnumeric(c));
    p.addParameter('CLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
    p.addParameter('Mask', [], @(m) isempty(m) || (islogical(m) && isequal(size(m), size(V))));
    p.addParameter('EdgeColor', 'none', @(s) ischar(s) || isstring(s) || (isnumeric(s) && numel(s) == 3));
    p.addParameter('XLabel', '', @(s) ischar(s) || isstring(s));
    p.addParameter('YLabel', '', @(s) ischar(s) || isstring(s));
    p.addParameter('Title',  '', @(s) ischar(s) || isstring(s));
    p.parse(nvArgs{:});
    opts = p.Results;

    % ----------------------
    % Construct edge coordinates from centers
    % ----------------------
    xedges = computeEdgesFromCenters(xcenters(:));
    yedges = computeEdgesFromCenters(ycenters(:).');

    % Build 2D edges grids
    Xedges2d = repmat(xedges(:), 1, numel(yedges));
    Yedges2d = repmat(yedges(:).', numel(xedges), 1);

    % Prepare CData with padded last row/col (not displayed by pcolor)
    C = nan(nx+1, ny+1);
    C(1:nx, 1:ny) = V;

    % Alpha mask
    if isempty(opts.Mask)
        alphaMask = ~isnan(C);
    else
        alphaMask = false(nx+1, ny+1);
        alphaMask(1:nx, 1:ny) = opts.Mask;
    end

    % ----------------------
    % Plot
    % ----------------------
    hpcolor = pcolor(ax, Xedges2d, Yedges2d, C);
    set(hpcolor, 'EdgeColor', opts.EdgeColor);
    set(hpcolor, 'AlphaData', double(alphaMask));

    % Colormap and CLim
    if ~isempty(opts.Colormap)
        colormap(ax, opts.Colormap);
    end
    if ~isempty(opts.CLim)
        caxis(ax, opts.CLim);
    end
    if opts.ShowColorbar
        colorbar(ax);
    end

    % Axis formatting
    if opts.AxisImage
        axis(ax, 'image');
    end
    xlabel(ax, string(opts.XLabel));
    ylabel(ax, string(opts.YLabel));
    title(ax,  string(opts.Title));

    % ----------------------
    % Centered ticks
    % ----------------------
    [xticks, xticklabels] = chooseTicks(xcenters(:), opts.NumTicks, opts.TickPrecision);
    [yticks, yticklabels] = chooseTicks(ycenters(:), opts.NumTicks, opts.TickPrecision);

    set(ax, 'XTick', xticks, 'XTickLabel', xticklabels);
    set(ax, 'YTick', yticks, 'YTickLabel', yticklabels);
end

% --------- Helpers ---------

function [xcenters, ycenters] = parseCentersXY(X, Y, nx, ny)
    % Accept vectors or meshgrid-like matrices
    if isvector(X)
        xcenters = reshape(X, [], 1);
    elseif ismatrix(X) && isequal(size(X), [nx, ny])
        % Meshgrid style: columns identical
        if ~all(all(abs(X - X(:,1)) < eps(max(abs(X(:))))))  % tolerance check
            error('pixplot:InvalidX', 'X must be a vector or meshgrid-like with identical columns.');
        end
        xcenters = X(:,1);
    else
        error('pixplot:InvalidX', 'Invalid X: must be vector length nx or matrix size [nx x ny].');
    end
    if numel(xcenters) ~= nx
        error('pixplot:InvalidXLength', 'Length of X (%d) must match size(V,1) (%d).', numel(xcenters), nx);
    end

    if isvector(Y)
        ycenters = reshape(Y, 1, []);
    elseif ismatrix(Y) && isequal(size(Y), [nx, ny])
        % Meshgrid style: rows identical
        if ~all(all(abs(Y - Y(1,:)) < eps(max(abs(Y(:))))))  % tolerance check
            error('pixplot:InvalidY', 'Y must be a vector or meshgrid-like with identical rows.');
        end
        ycenters = Y(1,:);
    else
        error('pixplot:InvalidY', 'Invalid Y: must be vector length ny or matrix size [nx x ny].');
    end
    if numel(ycenters) ~= ny
        error('pixplot:InvalidYLength', 'Length of Y (%d) must match size(V,2) (%d).', numel(ycenters), ny);
    end
end

function edges = computeEdgesFromCenters(centers)
    % Compute edges from center coordinates for possibly nonuniform spacing.
    n = numel(centers);
    if n == 1
        % Single pixel: infer unit spacing
        d = 1;
        edges = [centers(1) - 0.5*d; centers(1) + 0.5*d];
        return
    end
    d = diff(centers);
    edges = zeros(n+1,1);
    % Interior edges are midpoints
    edges(2:n) = centers(1:n-1) + d/2;
    % First and last edges extrapolated
    edges(1)   = centers(1)   - d(1)/2;
    edges(n+1) = centers(n)   + d(end)/2;
end

function [ticks, labels] = chooseTicks(centers, maxTicks, prec)
    n = numel(centers);
    step = max(1, round(n / maxTicks));
    idx = 1:step:n;
    ticks = centers(idx);
    fmt = sprintf('%%.%df', prec);
    labels = arrayfun(@(v) sprintf(fmt, v), ticks, 'UniformOutput', false);
end