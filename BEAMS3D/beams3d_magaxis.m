function [raxis, zaxis] = beams3d_magaxis(beam_data)
%BEASM3D_MAGAXIS Extracts the magnetic axis
%   The BEASM3D_MAGAXIS function returns R,Z position of the magnetic axis
%
%   Example:
%       beam_data=read_beams3d('beams3d_test.h5');
%       [raxis,zaxis] = beams3d_magaxis(beam_data);
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.0


raxis=[];
zaxis=[];
% Assume we want zeta=0 plane
phidex=1;

raxis=zeros(1,beam_data.nphi);
zaxis=zeros(1,beam_data.nphi);
% Make 2D
for i=1:beam_data.nphi-1
    S2D = squeeze(beam_data.S_ARR(:,i,:));
    smin = double(min(min(S2D)));
    [row,col] = find(S2D==smin);
    if all(diff(unique(col))==1) && all(diff(unique(row))==1)
        col=round(mean(col));
        row=round(mean(row));
    elseif numel(row)>1
        [c,~]=contour(beam_data.raxis,beam_data.zaxis,S2D',[1 1]);
        Z_min = [];
        idx = 1;
        N = size(c,2);
        [rg,zg]=ndgrid(beam_data.raxis,beam_data.zaxis);
        while idx < N
            pts = c(:,idx+(1:c(2,idx)));
            in_region = inpolygon(rg,zg,pts(1,:),pts(2,:));
            Z_min(end+1) = min(S2D(in_region));
            idx = idx+c(2,idx)+1;
        end
    end
    raxis(i) = beam_data.raxis(row);
    zaxis(i) = beam_data.zaxis(col); 
end
raxis(end) = raxis(1);
zaxis(end) = zaxis(1);

return;

end

