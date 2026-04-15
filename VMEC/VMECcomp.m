function output = VMECcomp( filename, datatype, varargin )
%VMECcomp(filename,datatype,plotmode) Creates plots comparing VMEC equilibria.
%   The VMECcomp function plots allows the user to create comparrision
%   plots of specific quantities between VMEC equilibria.  This is
%   accomplished by passing either a string with wildcards as the first
%   parameter (or a cell array of file names) and a string inicating the
%   type of value to plot as the second parameter.  The optional third
%   parameter selects the plot mode:
%       '3d'    3D waterfall plot (default, backwards compatible)
%       '2d'    Overlaid 2D line plot with legend
%   Available options are:
%       'curtor'        Toroidal Current
%       'extcur'        Vacuum Field Coil Currents
%       'iota'          Rotational Transform profile
%       'ac_aux'        Current spline
%       'am_aux'        Pressure spline
%       'q'             Safety Factor
%       'pressure'      Pressure profile
%       'current'       Toroidal current profile ('jcurv')
%       'jdotb'         <J*B>
%       'omega'         Rotation
%       'iota_press'    Pressure as a function of iota
%       'iota_pprime'   dp/ds as a function of iota
%       'flux'          Flux surfaces at phi=0
%       'fluxpi2'       Flux surfaces at quarter field period
%       'fluxpi'        Flux surfaces at half field period
%       'flux_edge'     VMEC edge and axis only phi=0
%       'magaxis'       3D plot of magnetic axis trajectory
%
% Example usage
%      haxis=VMECcomp('wout*','iota');       % 3D comparrision plot of iota
%      haxis=VMECcomp('wout*','iota','2d');  % 2D overlay plot of iota
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       2.0


use2d = false;
scale=1.0;
s_vals=[];
if nargin > 2
    i = 1;
    while i < numel(varargin)
        switch lower(varargin{i})
            case 'plotmode'
                i=i+1;
                use2d = strcmpi(varargin{i}, '2d');
            case 'scale'
                i=i+1;
                scale=varargin{i};
            case 's_vals'
                i=i+1;
                s_vals=varargin{i};
        end
        i=i+1;
    end
end



if ~iscell(filename)
    file_struct=dir(filename);
    nfiles=length(file_struct);
    filename=cell(1,nfiles);
    vmec_data=cell(1,nfiles);
    for i=1:nfiles
        filename{i}=file_struct(i).name;
        try
            vmec_data{i}=read_vmec(filename{i});
        catch
            vmec_data{i}=[];
        end
    end
else
    if isstr(filename{1})
        nfiles=max(size(filename));
        vmec_data=cell(1,nfiles);
        for i=1:nfiles
            try
                vmec_data{i}=read_vmec(filename{i});
            catch
                vmec_data{i}=[];
            end
        end
    elseif isstruct(filename{1})
        nfiles=max(size(filename));
        vmec_data=cell(1,nfiles);
        for i=1:nfiles
            vmec_data{i}=filename{i};
            filename{i} = strtrim(vmec_data{i}.input_extension);
        end
    end
    
end

if isempty(s_vals)
s_vals=[0.25, 0.5 0.75].^2;
end
s_label='s = \rho^2 [-]';
% % Fix underscores in filesname
% for i=1:length(filename)
%     filename{i} = strrep(filename{i},'_','\_');
% end

switch datatype
    case 'curtor'
        hold on
        curtor = [];
        for i=1:nfiles
            if ~isfield(vmec_data{i},'ctor'), continue; end
            curtor=[curtor; vmec_data{i}.ctor];
        end
        plot(1:nfiles,curtor,'o');
        set(gca,'XTick',1:nfiles);
        set(gca,'XTickLabel',filename);
        try; rotateXLabels(gca,90);end;
        ylabel('Net Toroidal Current');
        output=gca;
    case 'extcur'
        hold on
        extcur = [];
        for i=1:nfiles
            if ~isfield(vmec_data{i},'extcur'), continue; end
            ns=vmec_data{i}.ns;
            extcur=[extcur; vmec_data{i}.extcur];
        end
        if use2d
            plot(extcur');
            xlabel('Current Group');
            ylabel('Vaccum Field Currents');
            legend(filename,'Interpreter','none');
        else
            bar3(extcur);
            xlabel('Current Group');
            set(gca,'YTick',1:nfiles);
            set(gca,'YTickLabel',filename);
            zlabel('Vaccum Field Currents');
            view(3);
        end
        output=gca;
    case 'iota'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'iotaf'), continue; end
            ns=vmec_data{i}.ns;
            xdata{i} = 0:1/(ns-1):1;
            ydata{i} = vmec_data{i}.iotaf;
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            s_label, 'Rotational Transform \iota [-]', use2d);
        output=gca;
    case 'iota_flip'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'iotaf'), continue; end
            ns=vmec_data{i}.ns;
            xdata{i} = 0:1/(ns-1):1;
            ydata{i} = -vmec_data{i}.iotaf;
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            's = \rho^2', 'Rotational Transform \iota [-]', use2d);
        output=gca;
    case 'omega'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'omega'), continue; end
            ns=vmec_data{i}.ns;
            xdata{i} = 0:1/(ns-1):1;
            ydata{i} = vmec_data{i}.omega;
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            s_label, 'Rotational Transform', use2d);
        output=gca;
    case 'q'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'itoaf'), continue; end
            ns=vmec_data{i}.ns;
            xdata{i} = 0:1/(ns-1):1;
            ydata{i} = 1./vmec_data{i}.iotaf;
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            s_label, 'Safety Factor (q) [-]', use2d);
        output=gca;
    case 'pressure'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'presf'), continue; end
            ns=vmec_data{i}.ns;
            xdata{i} = 0:1/(ns-1):1;
            ydata{i} = vmec_data{i}.presf;
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            s_label, 'Pressure [Pa]', use2d);
        output=gca;
    case 'pressure_norm'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'presf'), continue; end
            ns=vmec_data{i}.ns;
            xdata{i} = 0:1/(ns-1):1;
            ydata{i} = vmec_data{i}.presf/vmec_data{i}.presf(1);
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            s_label, 'Normalized Pressure [-]', use2d);
        output=gca;        
    case 'jdotb'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'jdotb'), continue; end
            ns=vmec_data{i}.ns;
            xdata{i} = 0:1/(ns-1):1;
            ydata{i} = vmec_data{i}.jdotb;
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            s_label, '<J\cdotB>', use2d);
        output=gca;
    case 'buco'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'buco'), continue; end
            ns=vmec_data{i}.ns;
            xdata{i} = 0:1/(ns-1):1;
            ydata{i} = vmec_data{i}.buco;
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            s_label, '<B_u>', use2d);
        output=gca;
    case 'bvco'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'bvco'), continue; end
            ns=vmec_data{i}.ns;
            xdata{i} = 0:1/(ns-1):1;
            ydata{i} = vmec_data{i}.bvco;
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            s_label, '<B_V>', use2d);
        output=gca;
    case 'jcurv'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'jcurv'), continue; end
            ns=vmec_data{i}.ns;
            xdata{i} = 0:1/(ns-1):1;
            ydata{i} = vmec_data{i}.jcurv;
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            s_label, '<J^v>', use2d);
        output=gca;
    case 'jcurv_phi'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'jcurv'), continue; end
            xdata{i} = vmec_data{i}.phi;
            ydata{i} = vmec_data{i}.jcurv;
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            s_label, '<J^v>', use2d);
        output=gca;
    case 'jcuru'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'jcuru'), continue; end
            ns=vmec_data{i}.ns;
            xdata{i} = 0:1/(ns-1):1;
            ydata{i} = vmec_data{i}.jcuru;
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            s_label, '<J^u>', use2d);
        output=gca;
    case 'current'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'jcurv'), continue; end
            ns=vmec_data{i}.ns;
            xdata{i} = 0:1/(ns-1):1;
            ydata{i} = vmec_data{i}.jcurv.*2.*pi;
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            s_label, 'dI/ds [A]', use2d);
        output=gca;
    case 'iota_press'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'iotaf'), continue; end
            xdata{i} = vmec_data{i}.iotaf;
            ydata{i} = vmec_data{i}.presf;
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            'Rotational Transform [-]', 'Pressure [Pa]', use2d);
        output=gca;
    case 'iota_pprime'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'iotaf'), continue; end
            xdata{i} = vmec_data{i}.iotaf;
            ydata{i} = gradient(vmec_data{i}.presf, vmec_data{i}.phipf(1));
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            'Rotational Transform', 'dp/dpsi', use2d);
        output=gca;
    case 'ac_aux'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'acauxf'), continue; end
            ns = find(vmec_data{i}.acauxs > 0.0,1,'last');
            xdata{i} = vmec_data{i}.acauxs(1:ns);
            ydata{i} = vmec_data{i}.acauxf(1:ns);
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            s_label, 'Current Spline Coefficients', use2d);
        output=gca;
    case 'am_aux'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'amauxf'), continue; end
            ns = find(vmec_data{i}.amauxs > 0.0,1,'last');
            xdata{i} = vmec_data{i}.amauxs(1:ns);
            ydata{i} = vmec_data{i}.amauxf(1:ns);
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            s_label, 'Pressure Spline Coefficients', use2d);
        output=gca;
    case {'flux0','flux'}
        ntheta=90;
        hold on
        for i=1:nfiles
             plot(NaN, NaN,'Color',get_line_color(i, nfiles))
        end
        for i=1:nfiles
            if ~isfield(vmec_data{i},'rmnc'), continue; end
            color = get_line_color(i, nfiles);
            ns=vmec_data{i}.ns;
            rmnc=vmec_data{i}.rmnc;
            zmns=vmec_data{i}.zmns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            r=cfunct(0:2*pi/(ntheta-1):2*pi,0,rmnc,xm,xn);
            z=sfunct(0:2*pi/(ntheta-1):2*pi,0,zmns,xm,xn);
            if (vmec_data{i}.iasym)
                rmns=vmec_data{i}.rmns;
                zmnc=vmec_data{i}.zmnc;
                r=r+sfunct(0:2*pi/(ntheta-1):2*pi,0,rmns,xm,xn);
                z=z+cfunct(0:2*pi/(ntheta-1):2*pi,0,zmnc,xm,xn);
            end
            plot_flux_surfaces(r, z, ns, ntheta, i, s_vals, scale(i), color, use2d);
        end
        label_flux_axes(nfiles, filename, scale, use2d);
        axis equal
        output=gca;
    case {'flux_edge'}
        ntheta=90;
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'rmnc'), continue; end
            color = get_line_color(i, nfiles);
            ns=vmec_data{i}.ns;
            rmnc=vmec_data{i}.rmnc;
            zmns=vmec_data{i}.zmns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            r=cfunct(0:2*pi/(ntheta-1):2*pi,0,rmnc,xm,xn);
            z=sfunct(0:2*pi/(ntheta-1):2*pi,0,zmns,xm,xn);
            if (vmec_data{i}.iasym)
                rmns=vmec_data{i}.rmns;
                zmnc=vmec_data{i}.zmnc;
                r=r+sfunct(0:2*pi/(ntheta-1):2*pi,0,rmns,xm,xn);
                z=z+cfunct(0:2*pi/(ntheta-1):2*pi,0,zmnc,xm,xn);
            end
            if use2d
                plot(r(1,1),z(1,1),'+','Color',color);
                plot(r(ns,:),z(ns,:),'Color',color);
            else
                plot3(r(1,1),i,z(1,1),'+','Color',color);
                plot3(r(ns,:),i.*ones(1,ntheta),z(ns,:),color);
            end
        end
        label_flux_axes(nfiles, filename, scale, use2d);
        axis equal
        output=gca;
    case {'flux_edge3'}
        ntheta=90;
        zeta_vals = [0, 0, 0];
        % Pre-compute zeta values (need nfp from first valid file)
        for i=1:nfiles
            if isfield(vmec_data{i},'nfp')
                zeta_vals = [0, 2*pi/4/vmec_data{i}.nfp, pi/vmec_data{i}.nfp];
                break;
            end
        end
        for isub=1:3
            subplot(1,3,isub);
            hold on
            zeta = zeta_vals(isub);
            for i=1:nfiles
                if ~isfield(vmec_data{i},'rmnc'), continue; end
                color = get_line_color(i, nfiles);
                ns=vmec_data{i}.ns;
                rmnc=vmec_data{i}.rmnc;
                zmns=vmec_data{i}.zmns;
                xm=vmec_data{i}.xm;
                xn=vmec_data{i}.xn;
                r=cfunct(0:2*pi/(ntheta-1):2*pi,zeta,rmnc,xm,xn);
                z=sfunct(0:2*pi/(ntheta-1):2*pi,zeta,zmns,xm,xn);
                if (vmec_data{i}.iasym)
                    rmns=vmec_data{i}.rmns;
                    zmnc=vmec_data{i}.zmnc;
                    r=r+sfunct(0:2*pi/(ntheta-1):2*pi,zeta,rmns,xm,xn);
                    z=z+cfunct(0:2*pi/(ntheta-1):2*pi,zeta,zmnc,xm,xn);
                end
                if use2d
                    plot(r(1,1),z(1,1),'+','Color',color);
                    plot(r(ns,:),z(ns,:),'Color',color);
                else
                    plot3(r(1,1),i,z(1,1),'+','Color',color);
                    plot3(r(ns,:),i.*ones(1,ntheta),z(ns,:),color);
                end
            end
        label_flux_axes(nfiles, filename, scale, use2d);
            axis equal
        end
        output=gcf;
    case 'fluxpi2'
        ntheta=90;
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'rmnc'), continue; end
            color = get_line_color(i, nfiles);
            ns=vmec_data{i}.ns;
            rmnc=vmec_data{i}.rmnc;
            zmns=vmec_data{i}.zmns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            nfp=vmec_data{i}.nfp;
            r=cfunct(0:2*pi/(ntheta-1):2*pi,pi/2/nfp,rmnc,xm,xn);
            z=sfunct(0:2*pi/(ntheta-1):2*pi,pi/2/nfp,zmns,xm,xn);
            if (vmec_data{i}.iasym)
                rmns=vmec_data{i}.rmns;
                zmnc=vmec_data{i}.zmnc;
                r=r+sfunct(0:2*pi/(ntheta-1):2*pi,pi/2/nfp,rmns,xm,xn);
                z=z+cfunct(0:2*pi/(ntheta-1):2*pi,pi/2/nfp,zmnc,xm,xn);
            end
            plot_flux_surfaces(r, z, ns, ntheta, i, s_vals,scale(i), color, use2d);
        end
        label_flux_axes(nfiles, filename, scale, use2d);
        axis equal
        output=gca;
    case 'fluxpi'
        ntheta=90;
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'rmnc'), continue; end
            color = get_line_color(i, nfiles);
            ns=vmec_data{i}.ns;
            rmnc=vmec_data{i}.rmnc;
            zmns=vmec_data{i}.zmns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            nfp=vmec_data{i}.nfp;
            r=cfunct(0:2*pi/(ntheta-1):2*pi,pi/nfp,rmnc,xm,xn);
            z=sfunct(0:2*pi/(ntheta-1):2*pi,pi/nfp,zmns,xm,xn);
            if (vmec_data{i}.iasym)
                rmns=vmec_data{i}.rmns;
                zmnc=vmec_data{i}.zmnc;
                r=r+sfunct(0:2*pi/(ntheta-1):2*pi,pi/nfp,rmns,xm,xn);
                z=z+cfunct(0:2*pi/(ntheta-1):2*pi,pi/nfp,zmnc,xm,xn);
            end
            plot_flux_surfaces(r, z, ns, ntheta, i, s_vals,scale(i), color, use2d);
        end
        label_flux_axes(nfiles, filename, scale, use2d);
        axis equal
        output=gca;
    case 'magaxis'
        nzeta=360;
        cosph=cos(0:2*pi/(nzeta-1):2*pi);
        sinph=sin(0:2*pi/(nzeta-1):2*pi);
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'rmnc'), continue; end
            nfp=vmec_data{i}.nfp;
            ns=vmec_data{i}.ns;
            rmnc=vmec_data{i}.rmnc;
            zmns=vmec_data{i}.zmns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            r=cfunct(0,0:2*pi/(nzeta-1):2*pi,rmnc,xm,xn)*scale(i);
            z=sfunct(0,0:2*pi/(nzeta-1):2*pi,zmns,xm,xn)*scale(i);
            color = get_line_color(i, nfiles);
            if use2d
                plot(squeeze(r(1,1,:)).*cosph', ...
                     squeeze(r(1,1,:)).*sinph', 'Color', color);
            else
                plot3(squeeze(r(1,1,:)).*cosph', ...
                      squeeze(r(1,1,:)).*sinph', ...
                      squeeze(z(1,1,:)), color);
            end
        end
        xlabel('X [m]');
        ylabel('Y [m]');
        if use2d
            legend(filename,'Interpreter','none');
        else
            zlabel('Z [m]');
            view(3);
        end
        axis equal
        output=gca;
    case 'g'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'gmnc'), continue; end
            ns=vmec_data{i}.ns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            fmnc = vmec_data{i}.gmnc;
            f = cfunct(0,0,fmnc,xm,xn);
            if (vmec_data{i}.iasym)
                fmns = vmec_data{i}.gmns;
                f = f+sfunct(0,0,fmns,xm,xn);
            end
            xdata{i} = 0:1/(ns-1):1;
            ydata{i} = abs(f);
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            s_label, 'g', use2d);
        output=gca;
    case 'modb'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'bmnc'), continue; end
            ns=vmec_data{i}.ns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            fmnc = vmec_data{i}.bmnc;
            f = cfunct(0,0,fmnc,xm,xn);
            if (vmec_data{i}.iasym)
                fmns = vmec_data{i}.bmns;
                f = f+sfunct(0,0,fmns,xm,xn);
            end
            xdata{i} = 0:1/(ns-1):1;
            ydata{i} = abs(f);
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            s_label, '|B|', use2d);
        output=gca;
    case 'bsupu'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'bsupumnc'), continue; end
            ns=vmec_data{i}.ns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            fmnc = vmec_data{i}.bsupumnc;
            f = cfunct(0,0,fmnc,xm,xn);
            if (vmec_data{i}.iasym)
                fmns = vmec_data{i}.bsupumns;
                f = f+sfunct(0,0,fmns,xm,xn);
            end
            xdata{i} = 0:1/(ns-1):1;
            ydata{i} = abs(f);
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            s_label, 'B^U', use2d);
        output=gca;
    case 'bsupv'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'bsupvmnc'), continue; end
            ns=vmec_data{i}.ns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            fmnc = vmec_data{i}.bsupvmnc;
            f = cfunct(0,0,fmnc,xm,xn);
            if (vmec_data{i}.iasym)
                fmns = vmec_data{i}.bsupvmns;
                f = f+sfunct(0,0,fmns,xm,xn);
            end
            xdata{i} = 0:1/(ns-1):1;
            ydata{i} = abs(f);
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            s_label, 'B^V', use2d);
        output=gca;
    case 'special'
        hold on
        xdata = cell(1,nfiles); ydata = cell(1,nfiles);
        for i=1:nfiles
            if ~isfield(vmec_data{i},'nfp'), continue; end
            mu0=pi*4E-7;
            ns=vmec_data{i}.ns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            fmnc = vmec_data{i}.currvmnc;
            f = cfunct(0,0,fmnc,xm,xn);
            if (vmec_data{i}.iasym)
                fmns = vmec_data{i}.currvmns;
                f = f+sfunct(0,0,fmns,xm,xn);
            end
            xdata{i} = (0:1/(ns-1):1)';
            ydata{i} = abs(f).*mu0;
        end
        plot_profiles(xdata, ydata, nfiles, filename, ...
            s_label, 'currvmnc', use2d);
        output=gca;
    otherwise
        disp(['Error: Datatype ' strtrim(datatype) ' is not supported']);
        output = -1;
end

return

end

%--------------------------------------------------------------------------
% HELPER FUNCTIONS
%--------------------------------------------------------------------------

function plot_profiles(xdata, ydata, nfiles, filename, xlbl, zlbl, use2d)
%PLOT_PROFILES Plot profile data in 2D or 3D waterfall style.
%   Centralises the repeated plot3/plot logic used by most cases.
%   xdata, ydata  - cell arrays (1 x nfiles); empty cells are skipped.
%   filename       - cell array of display names (with escaped underscores).
%   xlbl, zlbl     - axis label strings (zlbl becomes ylabel in 2D mode).
%   use2d          - logical; true for overlaid 2D, false for 3D waterfall.

    if use2d
        for i=1:nfiles
            if isempty(xdata{i}), continue; end
            color = get_line_color(i, nfiles);
            lw = get_line_width(i, nfiles);
            plot(xdata{i}, ydata{i}, 'Color', color, 'LineWidth', lw);
        end
        xlabel(xlbl);
        ylabel(zlbl);
        legend(filename, 'Interpreter', 'none');
    else
        for i=1:nfiles
            if isempty(xdata{i}), continue; end
            ns = length(xdata{i});
            color = get_line_color(i, nfiles);
            lw = get_line_width(i, nfiles);
            plot3(xdata{i}, i.*ones(1,ns), ydata{i}, ...
                'Color', color, 'LineWidth', lw);
        end
        xlabel(xlbl);
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel(zlbl);
        view(3);
    end
end

function plot_flux_surfaces(r, z, ns, ntheta, idx, s_vals, scale, color, use2d)
%PLOT_FLUX_SURFACES Plot flux surface contours in 2D or 3D.
r=r*scale;
z=z*scale;
s_ind = round(s_vals * (ns - 1)) + 1;
plotst={'Color', color,'LineWidth',get_line_width(-1, 0)};
    if use2d
            plot(r(s_ind,:)', z(s_ind,:)',plotst{:});
        plot(r(1,1), z(1,1), '+',plotst{:});
        plot(r(ns,:), z(ns,:),plotst{:});
    else
        for j=s_ind
            plot3(r(j,:), idx.*ones(1,ntheta), z(j,:), color,'LineWidth',get_line_width(-1, 0));
        end
        plot3(r(1,1), idx, z(1,1), '+',plotst{:});
        plot3(r(ns,:), idx.*ones(1,ntheta), z(ns,:), color,'LineWidth',get_line_width(-1, 0));
    end
end

function label_flux_axes(nfiles, filename, scale,use2d)
%LABEL_FLUX_AXES Apply axis labels for flux surface plots.
    xlabel('R [m]');
    if any(scale~=1)
        % Append scale information to each filename entry for display
        % Ensure filename is a cell array of strings
        for k = 1:numel(filename)
            if scale(k)~=1
            filename{k} = sprintf('%s (scale=%.3g)', filename{k}, scale(k));
            end
        end
    end
    if use2d
        ylabel('Z [m]');
        legend(filename, 'Interpreter', 'none');
    else
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('Z [m]');
        view(3);
    end
end

function color = get_line_color(idx, nfiles)
%GET_LINE_COLOR Return line colour: blue for first, red for last, black otherwise.
    cols=lines(nfiles);
    color=cols(idx,:);
    % if idx == 1
    %     color = 'b';
    % elseif idx == nfiles
    %     color = 'r';
    % else
    %     color = 'k';
    % end
end

function lw = get_line_width(idx, nfiles)
%GET_LINE_WIDTH Return line width: thick for first/last, default otherwise.
    if idx == 1 || idx == nfiles
        lw = 0.5;% 2.0;
    else
        lw = 0.5;
    end
end
