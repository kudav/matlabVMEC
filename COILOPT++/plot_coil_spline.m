function plot_coil_spline( spline_data,varargin )
%PLOT_COIL_SPLINE(spline_data,options) Plots a spline file
% This function plots a spline file as generated by COILOPT++ and read by
% READ_COILOPT_SPLINE.
%
% Example usage
%      spline2coil('coil_spline_test.out','windsurf_test.nes');
%      coil_data=read_coils('coils.spline');
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.0

plottype='uv';
cwsfile = [];
lplotws = 0;
% Handle varargin
numdefargs=1;   %Number of default arguments
if nargin >numdefargs
    for i=1:nargin-numdefargs
        switch varargin{i}
            case {'uv','3d','3d_fp'}
                plottype=varargin{i};
            case {'cws'}
                cwsfile=varargin{i+1};
            case {'plot_cws'}
                lplotws = 1;
        end
    end
end

if ~isempty(cwsfile)
    fid = fopen(cwsfile,'r');
    if fid < 0
        disp(['Coil Winding Surface File Error: ' cwsfile]);
        return;
    end
    cwsdata=fscanf(fid,'%E',[4 inf]);
    fclose(fid);
    xm=cwsdata(1,:);
    xn=cwsdata(2,:);
    rbc=cwsdata(3,:);
    zbs=cwsdata(4,:);
    theta=0:2*pi/89:2*pi;
    zeta=0:2*pi/89:2*pi;
    r=squeeze(cfunct(theta,zeta,rbc',xm,xn));
    z=squeeze(sfunct(theta,zeta,zbs',xm,xn));
    if (lplotws)
        r_ws=cfunct(theta,zeta,rbc',xm,xn);
        z_ws=sfunct(theta,zeta,zbs',xm,xn);
    end
    u=theta/(2*pi);
    v=zeta/(2*pi);
    for i=1:length(zeta);
        x(:,i) = r(:,i).*cos(zeta(i)/spline_data.nfp);
        y(:,i) = r(:,i).*sin(zeta(i)/spline_data.nfp);
    end
end

coil_data.periods = spline_data.nfp;
datatype='coil_data';
coil_data.vert=[];
switch plottype
    case{'uv'}
        for i=1:spline_data.ncoil
            uvarr=fnplt(spline_data.sp{i});
        end
    case{'3d_fp'}
        for i=1:spline_data.ncoil
            coil_data.current_name{i} = ['MOD' num2str(i)];
            uvarr=fnplt(spline_data.sp{i});
            x_coil=interp2(u,v,x,uvarr(2,:),uvarr(1,:),'spline');
            y_coil=interp2(u,v,y,uvarr(2,:),uvarr(1,:),'spline');
            z_coil=interp2(u,v,z,uvarr(2,:),uvarr(1,:),'spline');
            verts=[x_coil; y_coil; z_coil];
            verts(4,:) = spline_data.current(i);
            verts(4,end) = 0.0;
            verts(5,:) = i;
            coil_data.vert = [coil_data.vert verts];
            coil_data.current_name{i} = ['MOD' num2str(i)];
            x_coil=interp2(u,v,x,1-uvarr(2,:),1-uvarr(1,:),'spline');
            y_coil=interp2(u,v,y,1-uvarr(2,:),1-uvarr(1,:),'spline');
            z_coil=interp2(u,v,z,1-uvarr(2,:),1-uvarr(1,:),'spline');
            verts=[x_coil; y_coil; z_coil];
            verts(4,:) = spline_data.current(i);
            verts(4,end) = 0.0;
            verts(5,:) = i;
            coil_data.vert = [coil_data.vert verts];
        end
    case{'3d'}
        for i=1:spline_data.ncoil
            coil_data.current_name{i} = ['MOD' num2str(i)];
            uvarr=fnplt(spline_data.sp{i});
            for np = 1:coil_data.periods
                x_coil=interp2(u,v,x,uvarr(2,:),uvarr(1,:),'spline');
                y_coil=interp2(u,v,y,uvarr(2,:),uvarr(1,:),'spline');
                z_coil=interp2(u,v,z,uvarr(2,:),uvarr(1,:),'spline');
                r_coil=sqrt(x_coil.*x_coil+y_coil.*y_coil);
                p_coil=atan2(y_coil,x_coil) + 2*pi*(np-1)./coil_data.periods;
                x_coil=r_coil.*cos(p_coil);
                y_coil=r_coil.*sin(p_coil);
                verts=[x_coil; y_coil; z_coil];
                verts(4,:) = spline_data.current(i);
                verts(4,end) = 0.0;
                verts(5,:) = i;
                coil_data.vert = [coil_data.vert verts];
                coil_data.current_name{i} = ['MOD' num2str(i)];
                x_coil=interp2(u,v,x,1-uvarr(2,:),1-uvarr(1,:),'spline');
                y_coil=interp2(u,v,y,1-uvarr(2,:),1-uvarr(1,:),'spline');
                z_coil=interp2(u,v,z,1-uvarr(2,:),1-uvarr(1,:),'spline');
                r_coil=sqrt(x_coil.*x_coil+y_coil.*y_coil);
                p_coil=atan2(y_coil,x_coil) + 2*pi*(np-1)./coil_data.periods;
                x_coil=r_coil.*cos(p_coil);
                y_coil=r_coil.*sin(p_coil);
                verts=[x_coil; y_coil; z_coil];
                verts(4,:) = spline_data.current(i);
                verts(4,end) = 0.0;
                verts(5,:) = i;
                coil_data.vert = [coil_data.vert verts];
            end
        end
end
switch plottype
    case{'uv'}
    case{'3d','3d_fp'}
        plot_coils(coil_data);
        if (lplotws)
            hold on;
            ha = isotoro(r_ws,z_ws,zeta./coil_data.periods,1);
            set(ha,'FaceAlpha',0.33,'FaceColor','blue');
        end

end
end

