function [P0_RphiZ, P1_RphiZ] = calc_beamvector(vmec_data, theta, zeta, alpha, beta, L, varargin)
% CALC_BEAMVECTOR Compute [R,phi,Z] points defining one or more beam lines from VMEC surface points.
%
% Usage:
%   [P0, P1] = calc_beamvector(vmec_data, theta, zeta, alpha, beta, L)
%   [P0, P1] = calc_beamvector(..., 'k', k_surface, 'degrees', true, 'plot', true)
%   [P0, P1] = calc_beamvector(..., 'addStellaratorSymmetry', true)
%
% Inputs (scalar or vectors):
%   vmec_data : structure as returned by READ_VMEC (must contain fields used below)
%   theta     : poloidal angle (VMEC u) [rad]. Scalar or vector
%   zeta      : toroidal angle (VMEC v; equals cylindrical phi) [rad]. Scalar or vector
%   alpha:    : rotation around the poloidal tangent e_pol. Positive alpha tilts the beam toward +e_tor.
%   beta:     : rotation around the toroidal tangent e_tor. Positive beta tilts the beam toward +e_pol.
%               The beam direction is obtained by starting from e_norm and applying the two rotations in the order alpha (about e_pol), then beta (about e_tor). Rotations are not commutative; if you prefer the opposite order, swap the two rotation blocks.
%   L         : distance between returned points, along the beam direction. Scalar or vector
%
% Name-value options:
%   'k'            : surface index (scalar or vector; 1..ns). Default: vmec_data.ns
%   'degrees'      : if true, alpha and beta are interpreted in degrees. Default: false (radians)
%   'plot'         : if true, plot the surface, points, and beam lines using isotoro. Default: false
%   'plotSurfaces' : vector of surface indices to show. Default: unique(k)
%   'ntheta'       : poloidal resolution for plotting. Default: 72
%   'nzeta'        : toroidal resolution for plotting. Default: 72
%   'axes'         : axes handle to plot into. Default: new figure
%   'surfaceColor' : RGB color for surface(s). Default: [0.7 0.7 0.8]
%   'lineColor'    : RGB color for beam line(s). Default: [1 0 0]
%   'p0Color'      : RGB color for P0 marker(s). Default: [0 0 0]
%   'p1Color'      : RGB color for P1 marker(s). Default: [0 0.6 0]
%   'lineWidth'    : line width for beam lines. Default: 2
%   'addStellaratorSymmetry' : if true, append beams at (theta, zeta + pi/nfp, alpha, beta) with
%                              theta -> -theta, alpha -> -alpha, beta -> -beta, L unchanged.
%
% Outputs:
%   P0_RphiZ  : N x 3 array of [R0, phi0, Z0] positions on the surface(s)
%   P1_RphiZ  : N x 3 array of [R1, phi1, Z1] positions L away along the beam direction
%
% Notes:
%   - Depends on cfunct, sfunct, and isotoro (for plotting).
%   - Beam direction: d = cos(alpha)*n + sin(alpha)*(cos(beta)*e1 + sin(beta)*e2),
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
addParameter(p, 'addStellaratorSymmetry', false);
parse(p, varargin{:});
k_in         = p.Results.k;
use_deg      = p.Results.degrees;
do_plot      = p.Results.plot;
plotS        = p.Results.plotSurfaces;
ntheta       = p.Results.ntheta;
nzeta        = p.Results.nzeta;
ax           = p.Results.axes;
surfaceColor = p.Results.surfaceColor;
lineColor    = p.Results.lineColor;
p0Color      = p.Results.p0Color;
p1Color      = p.Results.p1Color;
lineWidth    = p.Results.lineWidth;
addStelSym   = p.Results.addStellaratorSymmetry;

% Default to outermost surface if not given
if isscalar(k_in) && k_in == -1
    if isfield(vmec_data, 'ns')
        k_in = vmec_data.ns;
    else
        error('vmec_data.ns not found. Provide a valid k surface index.');
    end
end

% Ensure inputs are row vectors for broadcasting
theta = theta(:).';
zeta  = zeta(:).';
alpha = alpha(:).';
beta  = beta(:).';
L     = L(:).';

% Determine number of beams (before optional symmetry)
lens = [numel(theta), numel(zeta), numel(alpha), numel(beta), numel(L)];
N0 = max(lens);
if any(lens ~= 1 & lens ~= N0)
    error('theta, zeta, alpha, beta, and L must be scalars or vectors of the same length.');
end

% Expand scalars to length N0
theta = expand_to_len(theta, N0);
zeta  = expand_to_len(zeta,  N0);
alpha = expand_to_len(alpha, N0);
beta  = expand_to_len(beta,  N0);
L     = expand_to_len(L,     N0);

% k handling: scalar or vector of same length
if isscalar(k_in)
    k = repmat(k_in, 1, N0);
else
    k_in = k_in(:).';
    if numel(k_in) ~= N0
        error('If k is a vector, it must have the same length as the beam parameter arrays.');
    end
    k = k_in;
end

% Optionally add stellarator-symmetric beams (append to arrays)
if addStelSym
    if ~isfield(vmec_data, 'nfp') || isempty(vmec_data.nfp) || vmec_data.nfp <= 0
        error('addStellaratorSymmetry requested but vmec_data.nfp is missing or invalid.');
    end
    dphi_half = pi / vmec_data.nfp; % half a field period
    theta = [theta, -theta];
    zeta  = [zeta,  -zeta ];
    alpha = [alpha, -alpha];
    beta  = [beta,  -beta];
    L     = [L,      L];
    k     = [k,      k];
end

N = numel(theta);

% Degrees conversion for alpha/beta if needed
if use_deg
    alpha = alpha * pi/180;
    beta  = beta  * pi/180;
end

% Basic checks for k
if any(k < 1) || any(k > vmec_data.ns)
    error('Some requested surface indices k are out of range [1..%d].', vmec_data.ns);
end

% Prepare outputs
P0_RphiZ = nan(N, 3);
P1_RphiZ = nan(N, 3);

% Unpack VMEC Fourier data once
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
if iasym
    rmns  = vmec_data.rmns;
    zmnc  = vmec_data.zmnc;
    rumnc = vmec_data.rumnc;
    zumns = vmec_data.zumns;
    rvmnc = vmec_data.rvmnc;
    zvmns = vmec_data.zvmns;
end

% For plotting: store Cartesian endpoints
do_plot_beams = do_plot;
if do_plot_beams
    X0 = nan(N,1); Y0 = nan(N,1); Z0c = nan(N,1);
    X1 = nan(N,1); Y1 = nan(N,1); Z1c = nan(N,1);
end

% Compute each beam
for i = 1:N
    ti = theta(i);
    zi = zeta(i);
    ki = k(i);
    ai = alpha(i);
    bi = beta(i);
    Li = L(i);

    % Evaluate geometry and derivatives at (ti, zi) for all surfaces
    R_all = cfunct(ti, zi, rmnc,  xm, xn);
    Z_all = sfunct(ti, zi, zmns,  xm, xn);
    Ru_all = sfunct(ti, zi, rumns, xm, xn);
    Zu_all = cfunct(ti, zi, zumnc, xm, xn);
    Rv_all = sfunct(ti, zi, rvmns, xm, xn);
    Zv_all = cfunct(ti, zi, zvmnc, xm, xn);

    if iasym
        R_all  = R_all  + sfunct(ti, zi, rmns,  xm, xn);
        Z_all  = Z_all  + cfunct(ti, zi, zmnc,  xm, xn);
        Ru_all = Ru_all + cfunct(ti, zi, rumnc, xm, xn);
        Zu_all = Zu_all + sfunct(ti, zi, zumns, xm, xn);
        Rv_all = Rv_all + cfunct(ti, zi, rvmnc, xm, xn);
        Zv_all = Zv_all + sfunct(ti, zi, zvmns, xm, xn);
    end

    % Extract values at surface k
    R0 = squeeze(R_all(ki));
    Z0 = squeeze(Z_all(ki));
    Ru = squeeze(Ru_all(ki));
    Zu = squeeze(Zu_all(ki));
    Rv = squeeze(Rv_all(ki));
    Zv = squeeze(Zv_all(ki));

    phi0 = zi; % VMEC zeta equals cylindrical toroidal angle

    % Cartesian position of the base point
    x0 = R0*cos(phi0); y0 = R0*sin(phi0); z0 = Z0;

    % Tangent vectors (partials of X)
    Xt = [Ru*cos(phi0), Ru*sin(phi0), Zu];  % ∂X/∂theta
    Xz = [-R0*sin(phi0) + Rv*cos(phi0), ...
           R0*cos(phi0) + Rv*sin(phi0), ...
           Zv];                              % ∂X/∂zeta

   % Construct a local right-handed orthonormal frame:
    %   e_pol  = unit poloidal tangent (from Xt = ∂X/∂theta)
    %   e_tor  = unit toroidal tangent (orthogonalized from Xz = ∂X/∂zeta)
    %   e_norm = unit surface normal (e_pol × e_tor)
    tol = 1e-12;

    Xt_n = norm(Xt);
    Xz_n = norm(Xz);

    if Xt_n >= tol
        % Poloidal first, then toroidal via Gram-Schmidt
        e_pol = Xt / Xt_n;
        Xz_orth = Xz - dot(Xz, e_pol) * e_pol;
        Xz_orth_n = norm(Xz_orth);
        if Xz_orth_n < tol
            error('Toroidal tangent is degenerate at (theta,zeta) = (%.6g, %.6g) for beam %d.', ti, zi, i);
        end
        e_tor = Xz_orth / Xz_orth_n;
    elseif Xz_n >= tol
        % Toroidal first, then poloidal via Gram-Schmidt
        e_tor = Xz / Xz_n;
        Xt_orth = Xt - dot(Xt, e_tor) * e_tor;
        Xt_orth_n = norm(Xt_orth);
        if Xt_orth_n < tol
            error('Poloidal tangent is degenerate at (theta,zeta) = (%.6g, %.6g) for beam %d.', ti, zi, i);
        end
        e_pol = Xt_orth / Xt_orth_n;
    else
        error('Degenerate tangential basis at (theta,zeta) = (%.6g, %.6g) for beam %d.', ti, zi, i);
    end

    % Right-handed normal
    e_norm_vec = cross(e_pol, e_tor);
    e_norm_n = norm(e_norm_vec);
    if e_norm_n < tol
        error('Degenerate normal at (theta,zeta) = (%.6g, %.6g) for beam %d.', ti, zi, i);
    end
    e_norm = e_norm_vec / e_norm_n;

    % Beam direction from sequential rotations about e_pol and e_tor.
    % Start from the outward normal.
    v = e_norm;

    % Rotate around e_pol by alpha.
    % Convention: positive alpha tilts v toward +e_tor, hence the -sin(alpha) term.
    ca = cos(ai); sa = sin(ai);
    v = v*ca + cross(e_pol, v)*(-sa) + e_pol*dot(e_pol, v)*(1 - ca);

    % Then rotate around e_tor by beta.
    % Convention: positive beta tilts v toward +e_pol.
    cb = cos(bi); sb = sin(bi);
    v = v*cb + cross(e_tor, v)*(sb) + e_tor*dot(e_tor, v)*(1 - cb);

    % Normalize for numerical safety
    dn = norm(v);
    if ~isfinite(dn) || dn < tol
        error('Invalid beam direction at (theta,zeta) = (%.6g, %.6g) for beam %d.', ti, zi, i);
    end
    d = v / dn;

    % First point in Cartesian then to cylindrical
    x0 = x0 - 0.1*Li*d(1);
    y0 = y0 - 0.1*Li*d(2);
    z0 = z0 - 0.1*Li*d(3);    

    R0   = hypot(x0, y0);
    phi0 = atan2(y0, x0);
    Z0   = z0;        

    % Second point in Cartesian then to cylindrical
    x1 = x0 + Li*d(1);
    y1 = y0 + Li*d(2);
    z1 = z0 + Li*d(3);

    R1   = hypot(x1, y1);
    phi1 = atan2(y1, x1);
    Z1   = z1;



    % Pack outputs, normalize phi to [0, 2*pi)
    P0_RphiZ(i,:) = [R0, mod(phi0, 2*pi), Z0];
    P1_RphiZ(i,:) = [R1, mod(phi1, 2*pi), Z1];

    if do_plot_beams
        X0(i)=x0; Y0(i)=y0; Z0c(i)=z0;
        X1(i)=x1; Y1(i)=y1; Z1c(i)=z1;
    end
end

% Optional plotting using isotoro
if do_plot_beams
    if isempty(plotS)
        plotS = unique(k);
    end
    theta_vec = linspace(0, 2*pi, ntheta);
    zeta_vec  = linspace(0, 2*pi, nzeta);

    Rgrid = cfunct(theta_vec, zeta_vec, vmec_data.rmnc,  vmec_data.xm, vmec_data.xn);
    Zgrid = sfunct(theta_vec, zeta_vec, vmec_data.zmns,  vmec_data.xm, vmec_data.xn);
    if iasym
        Rgrid = Rgrid + sfunct(theta_vec, zeta_vec, vmec_data.rmns,  vmec_data.xm, vmec_data.xn);
        Zgrid = Zgrid + cfunct(theta_vec, zeta_vec, vmec_data.zmnc,  vmec_data.xm, vmec_data.xn);
    end

    if isempty(ax) || ~ishandle(ax)
        fig = figure('Color','w');
        ax = axes('Parent',fig);
    end
    axes(ax); hold(ax, 'on');

    try
        hsurf = isotoro(Rgrid, Zgrid, zeta_vec, plotS);
        try
            set(hsurf, 'FaceAlpha', 0.7, 'FaceColor', surfaceColor);
        catch
            for hh = reshape(hsurf,1,[])
                try, set(hh, 'FaceAlpha',0.7,'FaceColor',surfaceColor); end
            end
        end
    catch ME
        warning('isotoro plotting failed: %s. Proceeding with beam overlay only.', ME.message);
    end

    for i = 1:N
        plot3(ax, [X0(i) X1(i)], [Y0(i) Y1(i)], [Z0c(i) Z1c(i)], '-', 'Color', lineColor, 'LineWidth', lineWidth);
    end
    plot3(ax, X0, Y0, Z0c, 'o', 'MarkerSize', 8, 'MarkerFaceColor', p0Color, 'MarkerEdgeColor', 'k');
    plot3(ax, X1, Y1, Z1c, 's', 'MarkerSize', 8, 'MarkerFaceColor', p1Color, 'MarkerEdgeColor', 'k');

    axis(ax, 'equal');
    grid(ax, 'on');
    xlabel(ax, 'X');
    ylabel(ax, 'Y');
    zlabel(ax, 'Z');
    ttl = sprintf('Beam geometry on VMEC surface(s) k=%s', mat2str(unique(k)));
    title(ax, ttl);
    view(ax, 3);
end

end % main function

% Helper: expand scalars to target length
function v = expand_to_len(v, N)
    if isempty(v)
        v = zeros(1, N);
    elseif isscalar(v)
        v = repmat(v, 1, N);
    elseif numel(v) ~= N
        error('Input size mismatch.');
    end
end
