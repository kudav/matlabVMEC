function [P0_RphiZ, P1_RphiZ] = calc_beamvector(vmec_data, theta, zeta, alpha, beta, L, varargin)
% CALC_BEAMVECTOR Compute two [R,phi,Z] points defining a beam line from a VMEC surface point.
%
% Usage:
%   [P0, P1] = calc_beamvector(vmec_data, theta, zeta, alpha, beta, L)
%   [P0, P1] = calc_beamvector(..., 'k', k_surface, 'degrees', true, 'plot', true)
%
% Inputs:
%   vmec_data : structure as returned by READ_VMEC (must contain fields used below)
%   theta     : poloidal angle (VMEC u) at which to evaluate geometry (scalar)
%   zeta      : toroidal angle (VMEC v; equals cylindrical phi) (scalar)
%   alpha     : angle away from the local surface normal (0 = along normal)
%   beta      : rotation within tangent plane, measured from +∂X/∂theta direction
%   L         : distance between returned points, along the beam direction
%
% Name-value options:
%   'k'            : surface index (1..ns). Default: vmec_data.ns (outermost surface)
%   'degrees'      : if true, alpha and beta are interpreted in degrees. Default: false (radians)
%   'plot'         : if true, plot the surface, points, and the beam line using isotoro. Default: false
%   'plotSurfaces' : vector of surface indices to show. Default: k
%   'ntheta'       : poloidal resolution for plotting. Default: 72
%   'nzeta'        : toroidal resolution for plotting. Default: 72
%   'axes'         : axes handle to plot into. Default: new figure
%   'surfaceColor' : RGB color for surface(s). Default: [0.7 0.7 0.8]
%   'lineColor'    : RGB color for beam line. Default: [1 0 0]
%   'p0Color'      : RGB color for P0 marker. Default: [0 0 0]
%   'p1Color'      : RGB color for P1 marker. Default: [0 0.6 0]
%   'lineWidth'    : line width for beam line. Default: 2
%
% Outputs:
%   P0_RphiZ  : [R0, phi0, Z0] position on the surface at (theta, zeta)
%   P1_RphiZ  : [R1, phi1, Z1] position L away along the aimed beam direction
%
% Notes:
%   - This function depends on cfunct, sfunct, and isotoro (for plotting).
%   - The beam direction is: d = cos(alpha)*n + sin(alpha)*(cos(beta)*e1 + sin(beta)*e2),
%     where e1 is the normalized poloidal tangent (∂X/∂theta), and e2 = n × e1.

% Parse options
p = inputParser;
addParameter(p, 'k', -1);
addParameter(p, 'degrees', false);
addParameter(p, 'plot', false);
addParameter(p, 'plotSurfaces', []);
addParameter(p, 'ntheta', 72);
addParameter(p, 'nzeta', 72);
addParameter(p, 'axes', []);
addParameter(p, 'surfaceColor', [0.7 0.7 0.8]);
addParameter(p, 'lineColor', [1 0 0]);
addParameter(p, 'p0Color', [0 0 0]);
addParameter(p, 'p1Color', [0 0.6 0]);
addParameter(p, 'lineWidth', 2);
parse(p, varargin{:});
k           = p.Results.k;
use_deg     = p.Results.degrees;
do_plot     = p.Results.plot;
plotS       = p.Results.plotSurfaces;
ntheta      = p.Results.ntheta;
nzeta       = p.Results.nzeta;
ax          = p.Results.axes;
surfaceColor= p.Results.surfaceColor;
lineColor   = p.Results.lineColor;
p0Color     = p.Results.p0Color;
p1Color     = p.Results.p1Color;
lineWidth   = p.Results.lineWidth;

% Default to edge surface if not given
if k == -1
    if isfield(vmec_data, 'ns')
        k = vmec_data.ns;
        plotS=k;
    else
        error('vmec_data.ns not found. Provide a valid k surface index.');
    end
end

% Basic checks
if ~isscalar(theta) || ~isscalar(zeta) || ~isscalar(alpha) || ~isscalar(beta) || ~isscalar(L)
    error('theta, zeta, alpha, beta, and L must be scalars.');
end
if k < 1 || k > vmec_data.ns
    error('Requested surface index k=%d is out of range [1..%d].', k, vmec_data.ns);
end

% Unpack VMEC Fourier data
rmnc  = vmec_data.rmnc;
zmns  = vmec_data.zmns;
xm    = vmec_data.xm;
xn    = vmec_data.xn;
rumns = vmec_data.rumns;
zumnc = vmec_data.zumnc;
rvmns = vmec_data.rvmns;
zvmnc = vmec_data.zvmnc;

iasym = 0;
if isfield(vmec_data,'iasym'), iasym = vmec_data.iasym; end

% Evaluate geometry and derivatives at (theta, zeta), across all k, then pick index k
% R and Z at single point
R_all = cfunct(theta, zeta, rmnc,  xm, xn);
Z_all = sfunct(theta, zeta, zmns,  xm, xn);
% d/dtheta (u)
Ru_all = sfunct(theta, zeta, rumns, xm, xn);
Zu_all = cfunct(theta, zeta, zumnc, xm, xn);
% d/dzeta (v)
Rv_all = sfunct(theta, zeta, rvmns, xm, xn);
Zv_all = cfunct(theta, zeta, zvmnc, xm, xn);

if iasym
    % Asymmetric contributions
    rmns  = vmec_data.rmns;
    zmnc  = vmec_data.zmnc;
    rumnc = vmec_data.rumnc;
    zumns = vmec_data.zumns;
    rvmnc = vmec_data.rvmnc;
    zvmns = vmec_data.zvmns;

    R_all  = R_all  + sfunct(theta, zeta, rmns,  xm, xn);
    Z_all  = Z_all  + cfunct(theta, zeta, zmnc,  xm, xn);
    Ru_all = Ru_all + cfunct(theta, zeta, rumnc, xm, xn);
    Zu_all = Zu_all + sfunct(theta, zeta, zumns, xm, xn);
    Rv_all = Rv_all + cfunct(theta, zeta, rvmnc, xm, xn);
    Zv_all = Zv_all + sfunct(theta, zeta, zvmns, xm, xn);
end

% Extract values at the requested surface
R0 = squeeze(R_all(k));
Z0 = squeeze(Z_all(k));
Ru = squeeze(Ru_all(k));
Zu = squeeze(Zu_all(k));
Rv = squeeze(Rv_all(k));
Zv = squeeze(Zv_all(k));

phi0 = zeta; % VMEC zeta equals cylindrical toroidal angle in this convention

% Cartesian position of the base point
x0 = R0*cos(phi0); y0 = R0*sin(phi0); z0 = Z0;

% Tangent vectors (partials of X)
Xt = [Ru*cos(phi0), Ru*sin(phi0), Zu];  % ∂X/∂theta
Xz = [-R0*sin(phi0) + Rv*cos(phi0), ...
       R0*cos(phi0) + Rv*sin(phi0), ...
       Zv];                              % ∂X/∂zeta

% Normal and orthonormal frame
n = cross(Xt, Xz);
n_norm = norm(n);
if ~(isfinite(n_norm)) || n_norm < 1e-12
    error('Degenerate normal at (theta,zeta) = (%.6g, %.6g).', theta, zeta);
end
n_hat = n / n_norm;

t1 = Xt; t1n = norm(t1);
if t1n < 1e-12
    t1 = Xz; t1n = norm(t1);
end
if t1n < 1e-12
    error('Degenerate tangential basis at (theta,zeta) = (%.6g, %.6g).', theta, zeta);
end
e1 = t1 / t1n;
e2 = cross(n_hat, e1);
e2n = norm(e2);
if e2n < 1e-12
    error('Failed to construct orthonormal tangent frame.');
end
e2 = e2 / e2n;

% Angles: convert if needed
if use_deg
    alpha = alpha*pi/180;
    beta  = beta*pi/180;
end

% Beam direction in Cartesian
d = cos(alpha)*n_hat + sin(alpha)*(cos(beta)*e1 + sin(beta)*e2);
d = d / norm(d); % numeric safety

% Second point in Cartesian then to cylindrical
x1 = x0 + L*d(1);
y1 = y0 + L*d(2);
z1 = z0 + L*d(3);

R1   = hypot(x1, y1);
phi1 = atan2(y1, x1);
Z1   = z1;

% Pack outputs, normalize phi to [0, 2*pi) for convenience
P0_RphiZ = [R0, mod(phi0, 2*pi), Z0];
P1_RphiZ = [R1, mod(phi1, 2*pi), Z1];

% Optional plotting using isotoro
if do_plot
    % Choose surfaces to plot
    if isempty(plotS)
        plotS = k;
    end
    % Build theta and zeta grids for plotting
    theta_vec = linspace(0, 2*pi, ntheta);
    zeta_vec  = linspace(0, 2*pi, nzeta);

    % Compute r(s,theta,zeta) and z(s,theta,zeta) for all s on the grid
    Rgrid = cfunct(theta_vec, zeta_vec, rmnc,  xm, xn);
    Zgrid = sfunct(theta_vec, zeta_vec, zmns,  xm, xn);
    if iasym
        Rgrid = Rgrid + sfunct(theta_vec, zeta_vec, vmec_data.rmns,  xm, xn);
        Zgrid = Zgrid + cfunct(theta_vec, zeta_vec, vmec_data.zmnc,  xm, xn);
    end

    % Prepare axes
    if isempty(ax) || ~ishandle(ax)
        fig = figure('Color','w');
        ax = axes('Parent',fig);
    end
    axes(ax); hold(ax, 'on');

    % Plot equilibrium surface(s)
    try
        isotoro(Rgrid, Zgrid, zeta_vec, plotS, surfaceColor);
    catch ME
        warning('isotoro plotting failed: %s. Proceeding with beam overlay only.', ME.message);
    end

    % Overlay beam points and line
    plot3(ax, [x0 x1], [y0 y1], [z0 z1], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    plot3(ax, x0, y0, z0, 'o', 'MarkerSize', 8, 'MarkerFaceColor', p0Color, 'MarkerEdgeColor', 'k');
    plot3(ax, x1, y1, z1, 's', 'MarkerSize', 8, 'MarkerFaceColor', p1Color, 'MarkerEdgeColor', 'k');

    axis(ax, 'equal');
    grid(ax, 'on');
    xlabel(ax, 'X');
    ylabel(ax, 'Y');
    zlabel(ax, 'Z');
    title(ax, sprintf('Beam geometry on VMEC surface k=%d', k));
    view(ax, 3);
end
end