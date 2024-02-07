function [ lines_out ] = beams3d_apply_bpert(filename_in,fluxi0, ni,mi,varargin)
%beams3d_apply_bpert(filename_in,filename_out,fluxi0, ni,mi,varargin)
%applies a helical perturbation to the magnetic field in filename_in and
%writes the result to filename_out. Since only the magnetic field
%components are written in filename_out, this can also be a FIELDLINES
%restart file. The relative amplitude (to the maximum field found in the
%arrays) fluxi0, as well as the toroidal and poloidal mode numbers ni and
%mi can be specified as scalars or arrays, in which case all components are
%added up. Plotting of the perturbation amplitudes is available.
%
%The resulting magnetic field is calculated in the following way:
%The magnetic field in filename_in is scaled helically within the plasma
%according to fluxi0, ni and mi. The curl is calculated from the result,
%and this is taken as the perturbation field for B = B0 + B1.
%   Options:
%       'plot':        Plots the helical perturbation(s)
%
%   Usage:
%   vmec=read_vmec('wout_test.nc');
%   [pert] = beams3d_apply_bpert('beams3d_test.h5',1e-2,1,2,'plot','vmec',vmec);
%
%   Created by: D. Kulla (david.kulla@ipp.mpg.de)
%   Version:    1.0
%   Date:       22/11/30

lplot = 0;
lsave = 0;
llines=0;
lvmec=0;
lfieldlines=0;
calculation='strum'; %'hirv';
form='vacstrum';%'res','vacfcn'
phase=0;
% Handle varargin
if nargin > 4
    i = 1;
    while i <= numel(varargin)
        switch varargin{i}
            case 'plot'
                lplot=1;
            case 'fieldlines'
                lfieldlines=1;
                i=i+1;
                fieldlines_data = varargin{i};
            case 'save'
                lsave = 1;
                i=i+1;
                filename_out = varargin{i};
            case 'lines' %Calculate B_R and B_Z for fieldlines file
                llines=1;
            case 'vmec'
                lvmec = 1;
                i=i+1;
                vmec = varargin{i};
            case {'strum','hirv','ferrari'}
                calculation=varargin{i};
            case {'res','vacstrum','vacfcn'}
                form=varargin{i};       
            case {'phase'}
                i=i+1;
                phase = varargin{i};
        end
        i=i+1;
    end
end

br = h5read(filename_in,'/B_R');
bphi = h5read(filename_in,'/B_PHI');
bz = h5read(filename_in,'/B_Z');
r = h5read(filename_in,'/raxis');
phi = h5read(filename_in,'/phiaxis');
z = h5read(filename_in,'/zaxis');
rnorm=repmat(r,1,size(br,2),size(br,3));
% if llines
% br = br.*bphi./rnorm;
% bz = bz.*bphi./rnorm;
% end

sarr = h5read(filename_in,'/S_ARR');
uarr = h5read(filename_in,'/U_ARR');
sarr(sarr>1.5)=0;
u = zeros(size(sarr,1),size(sarr,3));
%uarr2 = zeros(size(sarr,1),size(sarr,3));
s=discretize(sarr,linspace(0,1,128));
s(isnan(s))=max(s,[],'all')+1;
if lvmec
    disp('USING PEST ANGLE FROM VMEC')
for i= 1:numel(r)
for j = 1:1
    for k = 1:numel(z)
        u(i,k)=vmec2pest(s(i,j,k),uarr(i,j,k),phi(j),vmec.lmns,vmec.xm,vmec.xn);
    end
end
end
uarr=permute(repmat(u,1,1,size(uarr,2)),[1,3,2]);
end

rhoarr = sqrt(sarr);
%[r0, phi0, z0] = meshgrid(r,phi,z);
[rg, phig, zg] = ndgrid(r,phi,z);
xg = rg .* cos(phig);
yg = rg .* sin(phig);
%bx = br .* cos(bphi);
%by = br.* sin(bphi);

switch form
    case 'res'
        % Given parameters
        alpha = 0.04;
        beta = 0.87;
        gamma = 0.01;
        s_21 = 0.2693;
        fluxir = islandWidthPerturbation(sarr, mi, ni, fluxi0, s_21, alpha, beta, gamma);
    case 'vacstrum'
        %Strumberger 2008, Perturbation 2: fluxi0=0.1
        fluxir= fluxi0.*(sarr).^(2/2) .* (1-sarr).^4;
        %Quadratic/analytical form in rho:
        %fluxir = fluxi0(i) * (rhoarr.^2) .*  (1- rhoarr).^2;
    case 'vacfcn'
        c1=2;
        c2=4;
        c3=1;
        f=@(x) x.^(c1/2).*(1-x).^c2;
        fluxir= fluxi0*f(sarr).^c3;
end
% 
% if lplot
%     figure
%     x=linspace(0,1,100);
%     plot(x,f(x));
%     hold on
%     %plot(x,fluxi0(i)*(sqrt(x)).^(2/2) .* (1-sqrt(x)).^4);
%     %plot(linspace(0,1,100),.1*linspace(0,1,100).^(2/2) .* (1-linspace(0,1,100)).^2);
%     %plot(rhoarr,fluxir);
%     %plot(rhoarr,fluxir_s);
%     xlabel('S=\rho^2')
%     ylabel('Perturbation Amplitude');
% end




%%
dr = r(2)-r(1);
dphi = phi(2)-phi(1);
dz = z(2)-z(1);

switch calculation
    case 'hirv'
        fluxphi = fluxir .* cos(mi.*uarr + ni .* phig-phase);
        fluxphi(sarr>1.5)=0;
        %Jacobsen formulation: B_pert = rot(alpha * B)
        brc = br.* fluxphi;
        bphic = bphi.* fluxphi;
        bzc = bz.* fluxphi;
        %Curl in cylindrical coordinates
        % Calculate the gradients of the components of B
        [~, dbrcdphi, dbrcdz] = gradient(brc, dr, dphi, dz);
        [~, ~, dbphicdz] = gradient(bphic, dr, dphi, dz);
        [dbzcdr, dbzcdphi, ~] = gradient(bzc, dr, dphi, dz);

        % Calculate the product r*bphic
        rbphic = rg .* bphic;

        % Now take the gradient of rbphic with respect to r
        [drbphicdr, ~, ~] = gradient(rbphic, dr, dphi, dz);

        % Calculate each component of the curl
        % Curl in cylindrical (r, phi, z) components
        curl_br = (1./rg).*dbzcdphi - dbphicdz;
        curl_bphi = dbrcdz - dbzcdr;
        curl_bz = (1./rg).*(drbphicdr - dbrcdphi);

        % Combine the components to form the curl vector
        bpert = cat(4, curl_br, curl_bphi, curl_bz);
    case 'strum'
        fluxphi = fluxir .* cos(mi.*uarr - ni .* phig-phase);
        fluxphi(sarr>1.5)=0;
        %Strumberger formulation B_pert = gradPsitilde x gradPhi
        %, dr, dphi, dz
        %[gradB, ~, ~] = calculateFieldProperties(lines_out);
        [dBdr, dBdphi, dBdz] = gradient(fluxphi, dr, dphi, dz);
        gradphi=permute(repmat(gradient(phi,dphi),1,numel(r),numel(z),3),[2,1,3,4]);
        gradB = cat(4, dBdr, dBdphi, dBdz); % 4th dimension represents the vector components
        bpert=cross(gradB,gradphi);
    case 'ferrari'
        fluxphi = fluxir .* cos(mi.*uarr + ni .* phig-phase);
        fluxphi(sarr>1.5)=0;        
        brc = br;%.* fluxphi;
        bphic = bphi.* fluxphi;
        bzc = bz;%.* fluxphi;        
        %[curl_br,curl_bphi,curl_bz,~] = curl(rg,phig,zg,brc,bphic,bzc);
        [curl_br,curl_bphi,curl_bz,~] =  curl(brc,bphic,bzc);
        bpert = cat(4, curl_br, curl_bphi, curl_bz);
end
curlr=squeeze(bpert(:,:,:,1));
curlphi=squeeze(bpert(:,:,:,2));
curlz=squeeze(bpert(:,:,:,3));
curlr(isnan(curlr)|isinf(curlr)|sarr>1) = 0;
curlphi(isnan(curlphi)|isinf(curlphi)|sarr>1) = 0;
curlz(isnan(curlz)|isinf(curlz)|sarr>1) = 0;

if lplot
    for j = numel(fluxi0)
        figure
        colors = parula(20);
        for i = 1:20
            p = patch(isosurface(xg,yg,zg,fluxphi,i/20*max(fluxphi(sarr<1),[],'all')));
            %isonormals(xg,yg,zg,fluxphi,p)
            p.FaceColor = colors(i,:);
            p.EdgeColor = 'none';
            hold on
        end
        view(3)
        axis equal
        camlight
        lighting phong
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
    end
    figure
    contour(r,z,squeeze(fluxphi(:,1,:))')
        xlabel('R')
        ylabel('Z')        
end

br = br + curlr;
bphi = bphi + curlphi;
bz = bz + curlz;

if llines
lines_out.B_R = br./bphi.*rnorm;
lines_out.B_Z = bz./bphi.*rnorm;
lines_out.datatype='FIELDLINES';
disp('Converted BR and BZ for FIELDLIENS!')
else
lines_out.B_R = br;
lines_out.B_Z = bz;
lines_out.datatype='OTHER';
end
lines_out.B_PHI = bphi;
lines_out.raxis = r;
lines_out.zaxis=z;
lines_out.phiaxis=phi;
lines_out.S_ARR=sarr;
lines_out.U_ARR=uarr;
lines_out.RHO_ARR=rhoarr;
lines_out.fluxphi=fluxphi;

lines_out.nr=numel(r);
lines_out.nphi=numel(phi);
lines_out.nz=numel(z);


if lsave
%end_state= h5read(filename_in,'/end_state');
%end_state=2.*ones(size(end_state));
%rbphi = h5read(filename_out,'/B_PHI');
%if sum(size(rbphi)-size(bphi))~=0
delete_hdf5_group(filename_out,'/B_R');
h5create(filename_out,'/B_R',size(br));
delete_hdf5_group(filename_out,'/B_PHI');
h5create(filename_out,'/B_PHI',size(br));
delete_hdf5_group(filename_out,'/B_Z');
h5create(filename_out,'/B_Z',size(br));
delete_hdf5_group(filename_out,'/raxis');
h5create(filename_out,'/raxis',size(r));
delete_hdf5_group(filename_out,'/phiaxis');
h5create(filename_out,'/phiaxis',size(phi));
delete_hdf5_group(filename_out,'/zaxis');
h5create(filename_out,'/zaxis',size(z));
delete_hdf5_group(filename_out,'/nr');
h5create(filename_out,'/nr',1);
delete_hdf5_group(filename_out,'/nphi');
h5create(filename_out,'/nphi',1);
delete_hdf5_group(filename_out,'/nz');
h5create(filename_out,'/nz',1);
%delete_hdf5_group(filename_out,'/end_state');
%h5create(filename_out,'/end_state',size(end_state));
%end
h5write(filename_out,'/B_R',br)
h5write(filename_out,'/B_PHI',bphi)
h5write(filename_out,'/B_Z',bz)
h5write(filename_out,'/raxis',r)
h5write(filename_out,'/phiaxis',phi)
h5write(filename_out,'/zaxis',z)
h5write(filename_out,'/nr',numel(r))
h5write(filename_out,'/nphi',numel(phi))
h5write(filename_out,'/nz',numel(z))
%h5write(filename_out,'/end_state',end_state)

end

if lfieldlines
[R,~,Z]=fieldlines_follow(lines_out,fieldlines_data.starts,fieldlines_data.phi_extent,fieldlines_data.poinc_loc,fieldlines_data.grid_extent);
figure
plot(R',Z','.')
end

end
% 
% function delete_hdf5_group(filename, group_name)
% fid = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
% H5L.delete(fid,group_name,'H5P_DEFAULT');
% H5F.close(fid);
% end
function A_mn = islandWidthPerturbation(s, m, n, rho_mn, s_mn, alpha, beta, gamma)
    % islandWidthPerturbation calculates the perturbation strength for a magnetic island.
    %
    % This function is vectorized to handle arrays of 's' values.
    %
    % Arguments:
    % s     : An array of normalized toroidal flux values of the q = m/n surface.
    % m     : Poloidal mode number.
    % n     : Toroidal mode number.
    % rho_mn: Free parameter for perturbation strength in the vacuum region.
    % s_mn  : Normalized toroidal flux at the q = m/n surface.
    % alpha : Free parameter alpha.
    % beta  : Free parameter beta.
    % gamma : Free parameter gamma.
    %
    % Returns:
    % A_mn  : An array of perturbation strengths for the magnetic island at each 's'.

    A_mn = zeros(size(s)); % Initialize A_mn with the same size as s

    % Logical array for indexing where s is less than or equal to s_mn
    index_s_less_equal = s <= s_mn;
    % Logical array for indexing where s is greater than s_mn
    index_s_greater = s > s_mn;

    % Calculate A_mn for s less than or equal to s_mn
    A_mn(index_s_less_equal) = rho_mn * alpha .* ((s(index_s_less_equal) / s_mn).^(m/2)) ...
                               .* (1 - beta * ((s(index_s_less_equal) / s_mn).^(1/2)));

    % Calculate A_mn for s greater than s_mn
    A_mn(index_s_greater) = rho_mn * (alpha * (1 - beta) ...
                               + gamma * (s(index_s_greater) / s_mn).^(1/2)) ...
                               ./ ((s(index_s_greater) / s_mn).^((n+1)/2));
end


