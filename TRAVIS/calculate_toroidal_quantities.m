function [W_MJ, beta_max, beta_vol] = calculate_toroidal_quantities(r_eff, Ne, Te_eV, R0, a, varargin)
%CALCULATE_TOROIDAL_QUANTITIES  Thermal energy and beta for a large‑aspect‑ratio torus.
%
%   INPUTS
%       r_eff   – radial coordinate (m), vector of length N, 0 <= r_eff <= a
%       Ne      – electron density profile (m^-3), same size as r_eff
%       Te_eV   – electron temperature profile (eV), same size as r_eff
%       R0      – major radius (m) (assumed constant for the whole plasma)
%       a       – minor radius (m) (upper limit of integration)
%       varargin – optional name/value pairs:
%           'B0'   – constant toroidal field (T).  If omitted, a field
%                    based on plasma current is used.
%           'Ip'   – plasma current (A).  Used only if B0 is not supplied.
%           'mu0'  – permeability of free space (default 4*pi*1e-7 H/m)
%
%   OUTPUTS
%       W_MJ       – total thermal energy (MJ)
%       beta_max   – maximum local beta (dimensionless)
%       beta_vol   – volume‑averaged beta (dimensionless)
%
%   NOTE
%       The function assumes a large‑aspect‑ratio tokamak so that
%           R(rho) ≈ R0
%       and the differential volume is
%           dV = 4*pi^2 * R0 * rho * drho .
%       If you need the full toroidal geometry, replace the dV definition
%       with the exact expression  dV = (R0 + rho*cos(theta))*rho d rho d theta d phi
%       and integrate over theta and phi analytically (gives the same 4*pi^2*R0 term).

% -------------------------------------------------------------------------
% 1.  Parse optional arguments
% -------------------------------------------------------------------------
p = inputParser;
addParameter(p,'B0',[]);          % constant toroidal field (T)
addParameter(p,'Ip',[]);          % plasma current (A)
addParameter(p,'mu0',4*pi*1e-7);   % H/m
parse(p,varargin{:});
B0   = p.Results.B0;
Ip   = p.Results.Ip;
mu0  = p.Results.mu0;

% -------------------------------------------------------------------------
% 2.  Basic checks
% -------------------------------------------------------------------------
assert(isvector(r_eff) && isvector(Ne) && isvector(Te_eV), ...
    'r_eff, Ne and Te_eV must be vectors of the same length.');
assert(numel(r_eff)==numel(Ne) && numel(Ne)==numel(Te_eV), ...
    'Input vectors must have identical length.');
assert(all(r_eff>=0) && all(r_eff<=a), 'r_eff must lie in [0 , a].');

% Ensure column vectors for broadcasting
r_eff = r_eff(:);
Ne    = Ne(:);
Te_eV = Te_eV(:);

% -------------------------------------------------------------------------
% 3.  Physical constants and unit conversion
% -------------------------------------------------------------------------
e   = 1.602176634e-19;   % elementary charge (C)
kB  = 1.380649e-23;      % Boltzmann constant (J/K)

% Convert temperature from eV to joules (T_J = T_eV * e)
Te_J = Te_eV * e;        % [J]

% -------------------------------------------------------------------------
% 4.  Magnetic field profile
% -------------------------------------------------------------------------
if ~isempty(B0)                     % user supplied constant field
    B = @(rho) B0 * ones(size(rho));
elseif ~isempty(Ip)                 % simple current‑driven field
    % B(phi) = mu0 * Ip / (2*pi*R)  (toroidal field from a circular loop)
    B = @(rho) mu0 * Ip ./ (2*pi*(R0 - rho));   % note: R ≈ R0 - rho for large aspect
else
    error('Either ''B0'' or ''Ip'' must be supplied.');
end

% -------------------------------------------------------------------------
% 5.  Differential volume element (large‑aspect‑ratio)
% -------------------------------------------------------------------------
dV = @(rho) 4*pi^2 * R0 .* rho;   % returns dV/drho (m^3/m)

% -------------------------------------------------------------------------
% 6.  Local pressure and beta
% -------------------------------------------------------------------------
p_local = Ne .* Te_J;               % [Pa] = N/m^2
beta_local = @(rho) (2*mu0 .* p_local ./ B(rho).^2);   % dimensionless

% -------------------------------------------------------------------------
% 7.  Thermal energy W = ∫ (3/2) p dV
% -------------------------------------------------------------------------
%integrand_W = @(rho) 1.5 .* p_local .* dV(rho);   % scalar for each rho
%W_J = integral(integrand_W, 0, a, 'ArrayValued', true);
%W_J = sum(1.5.*p_local.*dV(r_eff).*diff([0; r_eff]));
W_J = trapz(r_eff,1.5.*p_local.*dV(r_eff));
W_MJ = W_J * 1e-6;                               % convert to MJ

% -------------------------------------------------------------------------
% 8.  Maximum beta
% -------------------------------------------------------------------------
beta_vec = beta_local(r_eff);      % evaluate on the supplied grid
beta_max = max(beta_vec);

% -------------------------------------------------------------------------
% 9.  Volume‑averaged beta
% -------------------------------------------------------------------------
%integrand_beta = @(rho) beta_local(rho) .* dV(rho);
%beta_num = integral(integrand_beta, 0, a, 'ArrayValued', true);
%beta_den = integral(dV, 0, a, 'ArrayValued', true);
%beta_num = sum(beta_vec.*dV(r_eff).*diff([0; r_eff]));
beta_num = trapz(r_eff,beta_vec.*dV(r_eff));
%beta_den = sum(dV(r_eff).*diff([0; r_eff]));
beta_den = trapz(r_eff,dV(r_eff));
beta_vol = beta_num / beta_den;

% -------------------------------------------------------------------------
% 10.  Return results
% -------------------------------------------------------------------------
end