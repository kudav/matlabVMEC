function [r_out, v, z_out] = beams3d_suv2rzp(s, u, v, S_ARR, U_ARR, raxis, phiaxis, zaxis)
%BEAMS3D_SUV2RZP Converts (S, U, V) coordinates to (R, Z, Phi) coordinates.
%This is e.g. used in the volume normalization of the BEAMS3D distribution.
%   The function uses a Newton-Raphson method to iteratively solve for the 
%   (R, Z, Phi) coordinates from the given (S, U, V) coordinates, based on
%   provided grid arrays and axis information. Tolerance is set to 1e-9 and
%   the maximum number of iterations is 1000. If this is reached, a plot is
%   made on the result obtained and some debugging information is printed.
%
% Input Parameters:
%        s:            Normalized Toroidal Flux (S-coordinate)
%        u:            Poloidal Angle (U-coordinate)
%        v:            Toroidal Angle (V-coordinate)
%        S_ARR, U_ARR: Grid arrays for the S and U coordinates
%        raxis:        Major Radius (R) axis values
%        phiaxis:      Toroidal Angle (Phi) axis values
%        zaxis:        Vertical Distance (Z) axis values
%
% Output Parameters:
%        r_out:        Major Radius Distance R
%        phi_out:      Toroidal Angle Phi
%        z_out:        Vertical Distance Z
%
% Example Usage:
%        beam_data=read_beams3d('beams3d_test.h5');  
%        s = beam_data.dist_rhoaxis.^2;
%        u = beam_data.dist_uaxis;
%        v = beam_data.dist_paxis;
%        S_ARR = beam_data.S_ARR;
%        U_ARR = beam_data.U_ARR;
%        raxis = beam_data.raxis;
%        phiaxis = beam_data.phiaxis;
%        zaxis = beam_data.phiaxis;
%        [r,phi,z] = beams3d_suv2rzp(s,u,v,S_ARR,U_ARR,raxis,phiaxis,zaxis)
%
% Maintained by: David Kulla (david.kulla@ipp.mpg.de)
% Version:       1.00

% Create a grid from the input coordinates
[s, u, v] = ndgrid(sqrt(s), u, v);
[R_ARR, PHI_ARR, Z_ARR] = ndgrid(raxis, phiaxis, zaxis);

% Grid spacing for R, Phi, and Z axes
hr = raxis(2) - raxis(1);
hphi = phiaxis(2) - phiaxis(1);
hz = zaxis(2) - zaxis(1);

% Constants for the splines
tolerance = 1e-20;
% Constants and parameters for the Newton-Raphson method
nphi = numel(phiaxis);
pi2 = 2 * pi;
max_iterations = 1000;


% Initial values for R and Z using beams3d_magaxis function
data.nphi = nphi;
data.S_ARR = S_ARR;
data.raxis = raxis;
data.zaxis = zaxis;
[r_out, z_out] = beams3d_magaxis(data);

% Adjust PHI to be within the valid range
phi_out = mod(v, max(phiaxis));

% Interpolate R and Z values based on phi_out
r_out = interp1(phiaxis, r_out, squeeze(phi_out(1, 1, :)));
z_out = interp1(phiaxis, z_out, squeeze(phi_out(1, 1, :)));
r_out = permute(repmat(r_out, 1, size(phi_out, 1), size(phi_out, 2)), [2, 3, 1]);
z_out = permute(repmat(z_out, 1, size(phi_out, 1), size(phi_out, 2)), [2, 3, 1]);

% Adjust U to be within the valid range
u = mod(u, pi2);

% Convert S and U to Cartesian coordinates for the Newton-Raphson method
x0 = s .* cos(u);
y0 = s .* sin(u);

%fnorm = x0.^2 + y0.^2;
%fnorm = min(1 ./ fnorm, fnorm_max);
n = 1;
S_ARR=sqrt(S_ARR);
X_ARR=S_ARR.*cos(U_ARR);
Y_ARR=S_ARR.*sin(U_ARR);
interpmethod='linear';
X=griddedInterpolant(R_ARR,PHI_ARR,Z_ARR,X_ARR,interpmethod,'nearest');
Y=griddedInterpolant(R_ARR,PHI_ARR,Z_ARR,Y_ARR,interpmethod,'nearest');


[~, dXdR_ARR, dXdZ_ARR] = gradient(X_ARR,phiaxis,raxis,zaxis);
[~,dYdR_ARR, dYdZ_ARR] = gradient(Y_ARR,phiaxis,raxis,zaxis);

dXdR=griddedInterpolant(R_ARR,PHI_ARR,Z_ARR,dXdR_ARR,interpmethod,'nearest');
dXdZ=griddedInterpolant(R_ARR,PHI_ARR,Z_ARR,dXdZ_ARR,interpmethod,'nearest');
dYdR=griddedInterpolant(R_ARR,PHI_ARR,Z_ARR,dYdR_ARR,interpmethod,'nearest');
dYdZ=griddedInterpolant(R_ARR,PHI_ARR,Z_ARR,dYdZ_ARR,interpmethod,'nearest');
%figure
%hold on
% Loop (Newton's Method)
delR=zeros(size(r_out));
delZ=zeros(size(r_out));
residual=ones(size(r_out));
while ( max(residual,[],'all') > tolerance && n < max_iterations) %residual > tolerance &&

    x_term = x0 - X(r_out,phi_out,z_out);
    y_term = y0 - Y(r_out,phi_out,z_out);

    detJ = dXdR(r_out,phi_out,z_out) .* dYdZ(r_out,phi_out,z_out) - dYdR(r_out,phi_out,z_out) .* dXdZ(r_out,phi_out,z_out);
    detJ = max(abs(detJ), 0.01); % Upper bound for step size as detJ enters in denominator
    
    residual = (x_term.^2 + y_term.^2);%.* fnorm;
    dex=residual < 1e-4;

    delR = (x_term .* dYdZ(r_out,phi_out,z_out) - y_term .* dXdZ(r_out,phi_out,z_out)) ./ detJ;
    delZ = (-x_term .* dYdR(r_out,phi_out,z_out) + y_term .* dXdR(r_out,phi_out,z_out)) ./ detJ;

    delR = min(max(delR, -hr(1)), hr(1));
    delZ = min(max(delZ, -hz(1)), hz(1));

    delR(dex) = delR(dex) *0.5;% 1/sqrt(n);
    delZ(dex) = delZ(dex) * 0.5;%1/sqrt(n);

    r_out = max(min(r_out + delR, raxis(end)), raxis(1));
    z_out = max(min(z_out + delZ, zaxis(end)), zaxis(1));
    n = n + 1;

end

% If the maximum number of iterations is reached, print a warning.
if n >= max_iterations
    disp(['Residual max; ',num2str(max(residual,[],'all'))])
    figure
    scatter3(v(:),r_out(:),z_out(:),10.0,residual(:),'o');
    hold on
    [x_b3d,y_b3d,z_b3d] = meshgrid(phiaxis,raxis,zaxis);
    %tmp=permute(S_ARR,[3 1 2]);
    contourslice(x_b3d,y_b3d,z_b3d,S_ARR,[v(1,1,:)],[],[],linspace(0,1,10),interpmethod);%,linspace(-1,1,20)
    contourslice(x_b3d,y_b3d,z_b3d,U_ARR,[v(1,1,:)],[],[],linspace(0,2*pi,10),interpmethod);
    caxis([0 1])
    xlim([0 1.1])
    warning('Maximum number of iterations reached.');
end
end
