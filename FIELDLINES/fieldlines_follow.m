function [R_lines,PHI_lines,Z_lines] = fieldlines_follow(in_data,start_loc,phi_extent,poinc_loc,grid_extent)
%FIELDLINES_FOLLOW follows fieldlines from either a FIELDLINES or BEAMS3D
%run. It uses only the B_X components, so it can also be used e.g. from the
%beams3d_apply_bpert function
%
% Input parameters:
% in_data: data read by read_beams3d or read_fieldlines
% start_loc: num_lines by 2 array of R, Z starting positions
% phi_extent: toroidal range to integrate over
% poinc_loc: phi position to collect the poincare plot at
% grid_extent: R, Z grid extent
%
% Output parameters:
% R_lines: R coordinates at the poincare cross section
% PHI_lines: PHI coordinates at the poincare cross section
% Z_lines: Z coordinates at the poincare cross section
%
% Example usage:
% starts=[linspace(1.6,2.1,200)',repmat(0.03,200,1)];
% [pert] = beams3d_apply_bpert('beams3d_test.h5',5e-2,1,2,'ferrari','vacstrum');
% [R,PHI,Z]=fieldlines_follow(pert,starts,[0.0, 1000.0],2.0,[1.2,2.1,-.8,.8]);
%
% Maintained by: David Kulla (david.kulla@ipp.mpg.de)
% Version:       1.00
% Date  09/02/2024


nstart = size(start_loc,1);
[r,phi,z] = ndgrid(in_data.raxis,in_data.phiaxis,in_data.zaxis);
if numel(grid_extent)==0
    grid_extent(1)=in_data.raxis(1);
    grid_extent(2)=in_data.raxis(end);
    grid_extent(3)=in_data.zaxis(1);
    grid_extent(4)=in_data.raxis(end);
end
maxphi=max(in_data.phiaxis);
if ~strcmp(in_data.datatype,'FIELDLINES')
    dRdphi = r.*in_data.B_R ./ in_data.B_PHI;
    dZdphi = r.*in_data.B_Z ./ in_data.B_PHI;
else
    dRdphi=in_data.B_R;
    dZdphi=in_data.B_Z;
end


%Setup derivatives
dR_F = griddedInterpolant(r,phi,z,dRdphi,'cubic');
dZ_F = griddedInterpolant(r,phi,z,dZdphi,'cubic');


%Function to be integrated (q is a vector of length 2n with n entries for R
%and Z respectively.
    function dqdphi = f(phi,q)
        dqdphi = zeros(nstart*2,1);
        repphi=repmat(mod(phi,maxphi),nstart,1); %Periodic BC
        dqdphi(1:nstart) = dR_F(q(1:nstart),repphi,q(nstart+1:end));
        dqdphi(nstart+1:end)= dZ_F(q(1:nstart),repphi,q(nstart+1:end));
        dex=[q(1:nstart)<grid_extent(1) |...
            q(1:nstart)>grid_extent(2) ; ...
            q(nstart+1:end)<grid_extent(3) |...
            q(nstart+1:end)>grid_extent(4) ];
        dqdphi(dex)=0;
    end


    function [value,isterminal,direction] = events(t,~)
        value = mod(t-poinc_loc,maxphi)-pi;
        isterminal=0;
        direction=-1;
    end

options = odeset('RelTol',1e-5,'Events',@events);
starts = reshape(start_loc,1,[]);
[~,~,phie,ye,~] = ode89(@f,[phi_extent(1) phi_extent(2)],starts,options);
R_lines=ye(:,1:nstart)';
Z_lines=ye(:,nstart+1:end)';
PHI_lines=repmat(phie(:),1,nstart)';

end