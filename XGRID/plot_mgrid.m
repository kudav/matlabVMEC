function h = plot_mgrid(data, extcur, varargin)
% PLOT_MGRID(data, extcur, ...) Plot MGRID data with flexible options.
% 
% Compatible with the original usage:
%   plot_mgrid(data, extcur)                  % basic (single cutplane)
%   plot_mgrid(data, extcur, 'basic')
%   plot_mgrid(data, extcur, 'total')
%   plot_mgrid(data, extcur, '3dgrid')
%   plot_mgrid(data, extcur, 'modB')
%   plot_mgrid(data, extcur, 'cutplane', 3)   % cutplane index (1..nphi)
%
% New name-value options:
%   'PlotType'      : 'basic' | 'total' | '3dgrid' | 'modB'  (default 'basic')
%   'CutPlane'      : integer (phi index), default 1
%   'Figure'        : figure handle to draw in, default new figure
%   'Axes'          : axes handle to draw in (only used for single-axes plots)
%   'Colormap'      : colormap name or matrix, default 'parula'
%   'ColorLimits'   : [min max] for caxis, default auto per plot type
%   'ShowQuiver'    : logical, overlay quiver for br/bz (default true in basic/total)
%   'QuiverScale'   : numeric quiver scaling (default 1)
%   'Pause'         : seconds between frames for 'total' (default 0.5)
%   'GridResolution': integer resolution for 3d grid planes (default 20)
%
% Returns:
%   h : struct of handles (figure, axes, plots), fields depend on plot type.

% ---------------------------
% Input preprocessing
% ---------------------------
% Tolerate original positional flags and cutplane
[plotTypeFlag, cutplaneFlag, rest] = preprocessLegacyVarargin(varargin);

% Parse name-value options
p = inputParser;
p.addParameter('PlotType', plotTypeFlag, @(s)ischar(s) || isstring(s));
p.addParameter('CutPlane', cutplaneFlag, @(x)isscalar(x) && isnumeric(x) && x>=1);
p.addParameter('Figure', [], @(h) isempty(h) || ishghandle(h,'figure'));
p.addParameter('Axes', [], @(h) isempty(h) || ishghandle(h,'axes'));
p.addParameter('Colormap', 'parula', @(c) (ischar(c) || isstring(c) || isnumeric(c)));
p.addParameter('ColorLimits', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
p.addParameter('ShowQuiver', true, @(b)islogical(b) || isnumeric(b));
p.addParameter('QuiverScale', 1, @(x)isscalar(x) && isnumeric(x));
p.addParameter('Pause', 0.5, @(x)isscalar(x) && isnumeric(x) && x>=0);
p.addParameter('GridResolution', 20, @(x)isscalar(x) && isnumeric(x) && x>=2);
p.parse(rest{:});
opts = p.Results;

% ---------------------------
% Validate data and extcur
% ---------------------------
% Determine number of currents
if isfield(data,'nextcur')
    ncur = data.nextcur;
else
    % Infer from br size if possible
    if isfield(data,'br') && ndims(data.br)==4
        ncur = size(data.br,4);
    else
        error('Cannot determine number of external currents (nextcur).');
    end
end

% Accept scalar extcur (replicate), or vector matching ncur
if isempty(extcur)
    error('extcur cannot be empty.');
end
if isscalar(extcur) && ncur>1
    extcur = repmat(extcur, 1, ncur);
elseif numel(extcur) ~= ncur
    error('Extcur size mismatch: expected %d elements, got %d.', ncur, numel(extcur));
end

% Validate presence of required fields
reqFields = {'br','bphi','bz','raxis','zaxis'};
for k = 1:numel(reqFields)
    if ~isfield(data, reqFields{k})
        error('Data missing required field: %s', reqFields{k});
    end
end
% Determine nphi and phi vector
if isfield(data,'phi')
    phi = data.phi(:).';
else
    % Fall back to uniform spacing if missing
    nphi = size(data.br,3);
    phi = linspace(0, 2*pi, nphi);
end

% Dimensions
nr   = size(data.br,1);
nz   = size(data.br,2);
nphi = size(data.br,3);

% Validate cutplane index
cutIdx = max(1, min(nphi, round(opts.CutPlane)));

% ---------------------------
% Compute total fields efficiently
% ---------------------------
scale = reshape(extcur(:), [1 1 1 ncur]);  % broadcast along 4th dimension
brt   = sum(data.br   .* scale, 4);
bphit = sum(data.bphi .* scale, 4);
bzt   = sum(data.bz   .* scale, 4);
bmag  = sqrt(brt.^2 + bphit.^2 + bzt.^2);

% Axes grids consistent with nr x nz layout
[R, Z] = ndgrid(data.raxis(:), data.zaxis(:));  % size [nr x nz]

% Color limits defaults
autoClimComponents = [min([min(brt,[],'all'), min(bzt,[],'all'), min(bphit,[],'all')]), ...
                      max([max(brt,[],'all'), max(bzt,[],'all'), max(bphit,[],'all')])];
autoClimBmag = [min(bmag,[],'all'), max(bmag,[],'all')];

% ---------------------------
% Figure / axes setup
% ---------------------------
h = struct();
if ~isempty(opts.Axes)
    axParent = opts.Axes;
    h.figure = ancestor(axParent,'figure');
else
    if isempty(opts.Figure) || ~ishghandle(opts.Figure,'figure')
        h.figure = figure('Position',[100 100 1280 720]); % more reasonable default
    else
        h.figure = opts.Figure;
    end
    axParent = [];
end
colormap(h.figure, opts.Colormap);

% ---------------------------
% Plotting by type
% ---------------------------
plotType = lower(string(opts.PlotType));
switch plotType
    case "total"
        % Pan through all phi cuts, updating plots
        if isempty(axParent)
            ax1 = subplot(2,2,1,'Parent',h.figure);
            ax2 = subplot(2,2,2,'Parent',h.figure);
            ax3 = subplot(2,2,3,'Parent',h.figure);
            ax4 = subplot(2,2,4,'Parent',h.figure);
        else
            % If a single axes was provided, we will use it only for |B| plot
            ax1 = axParent; ax2 = axParent; ax3 = axParent; ax4 = axParent;
        end
        clim = chooseClim(opts.ColorLimits, autoClimComponents);
        for i = 1:nphi
            % Toroidal field
            axes(ax1);
            hp1 = pcolor(R, Z, bphit(:,:,i)); set(hp1,'EdgeColor','none');
            xlabel('R'); ylabel('Z'); title(sprintf('Toroidal B_\\phi, \\phi = %.1f°', rad2deg(phi(i))));
            colorbar; caxis(clim); axis image;

            % Radial field
            axes(ax2);
            hp2 = pcolor(R, Z, brt(:,:,i)); set(hp2,'EdgeColor','none');
            xlabel('R'); ylabel('Z'); title(sprintf('Radial B_r, \\phi = %.1f°', rad2deg(phi(i))));
            colorbar; caxis(clim); axis image;

            % Vertical field
            axes(ax3);
            hp3 = pcolor(R, Z, bzt(:,:,i)); set(hp3,'EdgeColor','none');
            xlabel('R'); ylabel('Z'); title(sprintf('Vertical B_z, \\phi = %.1f°', rad2deg(phi(i))));
            colorbar; caxis(clim); axis image;

            % Combined with quiver
            axes(ax4);
            hp4 = pcolor(R, Z, bphit(:,:,i)); set(hp4,'EdgeColor','none'); hold on;
            if opts.ShowQuiver
                quiver(R, Z, opts.QuiverScale*brt(:,:,i), opts.QuiverScale*bzt(:,:,i), 'k');
            end
            hold off; colorbar; caxis(clim);
            xlabel('R'); ylabel('Z'); title(sprintf('B-field components, \\phi = %.1f°', rad2deg(phi(i))));
            axis image;

            pause(opts.Pause);
        end
        h.axes = [ax1, ax2, ax3, ax4];

    case "basic"
        % Single cutplane (default)
        if isempty(axParent)
            ax1 = subplot(2,2,1,'Parent',h.figure);
            ax2 = subplot(2,2,2,'Parent',h.figure);
            ax3 = subplot(2,2,3,'Parent',h.figure);
            ax4 = subplot(2,2,4,'Parent',h.figure);
        else
            % Respect provided axes by drawing only the combined plot
            ax1 = axParent; ax2 = axParent; ax3 = axParent; ax4 = axParent;
        end
        clim = chooseClim(opts.ColorLimits, autoClimComponents);

        axes(ax1);
        hp1 = pcolor(R, Z, bphit(:,:,cutIdx)); set(hp1,'EdgeColor','none');
        xlabel('R'); ylabel('Z'); title('Toroidal B_\phi'); colorbar; caxis(clim); axis image;

        axes(ax2);
        hp2 = pcolor(R, Z, brt(:,:,cutIdx)); set(hp2,'EdgeColor','none');
        xlabel('R'); ylabel('Z'); title('Radial B_r'); colorbar; caxis(clim); axis image;

        axes(ax3);
        hp3 = pcolor(R, Z, bzt(:,:,cutIdx)); set(hp3,'EdgeColor','none');
        xlabel('R'); ylabel('Z'); title('Vertical B_z'); colorbar; caxis(clim); axis image;

        axes(ax4);
        hp4 = pcolor(R, Z, bphit(:,:,cutIdx)); set(hp4,'EdgeColor','none'); hold on;
        if opts.ShowQuiver
            quiver(R, Z, opts.QuiverScale*brt(:,:,cutIdx), opts.QuiverScale*bzt(:,:,cutIdx), 'k');
        end
        hold off; colorbar; caxis(clim);
        xlabel('R'); ylabel('Z'); title(sprintf('B components at \\phi = %.1f°', rad2deg(phi(cutIdx))));
        axis image;

        h.axes = [ax1, ax2, ax3, ax4];

    case "modb"
        % Plot |B| at selected cutplane
        ax = axParent;
        if isempty(ax)
            ax = axes('Parent', h.figure);
        end
        clim = chooseClim(opts.ColorLimits, autoClimBmag);
        hp = pcolor(R, Z, bmag(:,:,cutIdx)); set(hp,'EdgeColor','none');
        xlabel('R'); ylabel('Z'); title('|B|'); colorbar; caxis(clim); axis image;
        h.axes = ax; h.plots.modB = hp;

    case "3dgrid"
        % Visualize cutplanes as surfaces in 3D
        ax = axParent;
        if isempty(ax)
            ax = axes('Parent', h.figure);
        end
        hold(ax,'on');
        % Build a coarse grid for visualization
        rmin = pickField(data,'rmin', min(data.raxis));
        rmax = pickField(data,'rmax', max(data.raxis));
        zmin = pickField(data,'zmin', min(data.zaxis));
        zmax = pickField(data,'zmax', max(data.zaxis));
        rtemp = linspace(rmin, rmax, opts.GridResolution);
        ztemp = linspace(zmin, zmax, opts.GridResolution);
        [Rcoarse, Zcoarse] = ndgrid(rtemp, ztemp);
        for j = 1:nphi
            X = Rcoarse .* cos(phi(j));
            Y = Rcoarse .* sin(phi(j));
            Zgrid = Zcoarse;
            s = surf(ax, X, Y, Zgrid, 'FaceColor', 'none', 'EdgeColor', [0.5 0.5 0.5]);
        end
        hold(ax,'off');
        xlabel(ax,'X'); ylabel(ax,'Y'); zlabel(ax,'Z');
        axis(ax,'equal'); xlim(ax,[0 rmax]); ylim(ax,[-rmax rmax]); zlim(ax,[zmin zmax]);
        view(ax,3);
        h.axes = ax;

    otherwise
        error('Unknown PlotType: %s', opts.PlotType);
end

end

% ---------------------------
% Helper functions
% ---------------------------
function [plotType, cutplane, rest] = preprocessLegacyVarargin(args)
% Accept legacy positional flags: 'basic'|'total'|'3dgrid'|'modB' and 'cutplane', val
plotType = 'basic';
cutplane = 1;
rest = args;
if isempty(args), return; end
% Find plot type flag
flags = {'basic','total','3dgrid','modB'};
for k = 1:numel(args)
    if ischar(args{k}) || isstring(args{k})
        s = lower(string(args{k}));
        if any(strcmpi(s, flags))
            plotType = char(s);
            rest(k) = []; % remove flag
            break;
        end
    end
end
% Find 'cutplane' followed by a value
idx = [];
for k = 1:numel(rest)
    if ischar(rest{k}) || isstring(rest{k})
        if strcmpi(string(rest{k}), 'cutplane') && (k < numel(rest))
            val = rest{k+1};
            if isnumeric(val) && isscalar(val)
                cutplane = val;
                idx = [k k+1];
            end
            break;
        end
    end
end
if ~isempty(idx)
    rest(idx) = []; % remove cutplane pair
end
end

function clim = chooseClim(userClim, autoClim)
if isempty(userClim), clim = autoClim; else, clim = userClim; end
end

function val = pickField(s, name, default)
if isfield(s, name), val = s.(name); else, val = default; end
end
