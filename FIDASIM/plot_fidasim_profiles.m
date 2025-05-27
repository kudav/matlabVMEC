function plot_data = plot_fidasim_profiles(filename,in_data,varargin)
%PLOT_FIDASIM_PROFILES plots radial profiles from the FIDASIM code. As
%these are used to compare to experimental data, the integration range and
%dispersion is set by variables in the in_data structure generated e.g. by
%the geto_bes_fida_aug_data function. Plots are shared between the
%two functions.
%
% Example usage
%      [~,in_data] = get_bes_fida_aug_data(filename,in_data,'fidabes');
%      plot_fidasim_profiles(filename,in_data,'X');
%      !!! 'X' can be 'fida', 'bes', or 'fidabes'
%      plot_fidasim_profiles(filename,_,'spec_bes'); %Forces calculating BES from
%      total spectrum, not from full energy component
%
% Miscellaneous Arguments
%      plot_fidasim(runid,'mean'); %Apply moving mean to spectrum
%      plot_fidasim(runid,'in_data',in_data); %Used for dispersion
%      plot_fidasim(runid,'save'); %Export figures (.fig and .png)
%      plot_fidasim(runid,'name', 'test'); %ID Name for legend
%      plot_fidasim(runid,'fac', 1.0); %Scaling factor
%
if isfield(in_data,'bg_range')
    bg_range=in_data.bg_range;
else
    bg_range = [664.5, 666];
end
if isfield(in_data,'bes_range')
    bes_range = in_data.bes_range;
else
    bes_range=[];
end
if isfield(in_data,'fida_range')
    fida_range = in_data.fida_range;
else
    fida_range = [659.5, 660.5];
end



lsave = 0;
lmean = 0;
lload_fidasim=0;
lrho=0;
liota=0;%to plot iota values
efit={};
vmec={};
ax={};
leps=0;
lxerr=0;%plot x error from neutral density
lspecbes=0; %calculate BES from full spectrum
lrunidfromdat=0;%to use runid from inputs dat for results
plot_type = {};
data={};%FIDASIM results data (spec)
channel_spec=0;
geom={};
b3d_name='';
linestyle = '+';
fac = 1;

if strcmp(filename(end-10:end),'_spectra.h5')
    filename=filename(1:end-11);
end
name = filename(1:end-2);
filename_runid=filename;

if nargin > 2
    i = 1;
    while i < nargin-1
        switch varargin{i}
            case {'FIDA','BES','FIDABES','fida','bes','fidabes','bck',...
                    'fidaspec','fida_bck','fidabck','fidabg','bes_comp'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
            case 'spectrum'
                plot_type{end+1}=varargin{i};
                i=i+1;
                if ischar(varargin{i}) %TODO: Allow for multiple channels to be specified
                    channel_spec = find(strcmp(varargin{i},cellstr(deblank(names'))));
                elseif iscell(varargin{i})
                    for j = 1:numel(varargin{i})
                        channel_spec(j) = find(strcmp(varargin{i}{j},cellstr(deblank(names'))));
                    end
                else
                    channel_spec = varargin{i};
                end
            case 'mean'
                lmean =1;
            case 'save'
                lsave = 1;
            case 'avg_frames'
                i = i+1;
                avg_frames = varargin{i};
            case 'spec_bes'
                lspecbes=1;
            case 'rho'
                lrho=1;
                lgeom=1;
                i=i+1;
                b3d_name=varargin{i};
            case 'eqdsk'
                lrho=1;
                lgeom=1;
                i=i+1;
                efit=read_efit(varargin{i});
            case 'vmec'
                lrho=1;
                lgeom=1;
                if ischar(varargin{i})&~isempty(vmec)
                    i=i+1;
                    vmec=read_vmec(['wout_',varargin{i},'.nc']);
                else
                    i=i+1;
                    vmec=varargin{i};
                end
            case 'qeqdsk'
                liota=1;
                lgeom=1;
                i=i+1;
                efit=read_efit(varargin{i});
                i=i+1;
                iota_vec=varargin{i};
            case 'qvmec'
                liota=1;
                i=i+1;
                if ischar(varargin{i})&~isempty(vmec)
                    vmec=read_vmec(['wout_',varargin{i},'.nc']);
                    i=i+1;
                else
                    vmec=varargin{i};
                    i=i+1;
                end
                iota_vec=varargin{i};
            case {'frominputs','fromdat'}
                lrunidfromdat=1;
            case 'ax'
                i=i+1;
                ax = varargin{i};
            case 'fac'
                i = i+1;
                fac = varargin{i};
            case 'load'
                lload_fidasim=1;
                i=i+1;
                channel= varargin{i};
            case 'xerr'
                lxerr=1;
                lrunidfromdat=1;
            case 'name'
                i = i+1;
                name = varargin{i};
            case 'eps'
                leps=1;
            case 'style'
                i = i+1;
                linestyle = varargin{i};
            otherwise
                disp(['ERROR: Option ', varargin{i}, ' not found!']);

        end
        i=i+1;
    end
end

if isfield(in_data,'dex')
    dex = in_data.dex;
else
    dex = [];
end


if isfield(in_data,'ax') &&isempty(ax)
    ax=in_data.ax;
end
if lrunidfromdat
    nml_name=[filename,'_inputs.dat'];
    if isfile(nml_name)
        input=read_namelist(nml_name,'fidasim_inputs');
        filename_runid=input.runid;
    else
        disp(['Could not find namelist with filename ', nml_name])
        return
    end
end

spec_name = [filename_runid,'_spectra.h5'];
geom_name = [filename,'_geometry.h5'];

%if ~isfile(geom_name)|| (lload_fidasim && ~isfile(spec_name))||lrho
    if lrho
        if numel(b3d_name)~=0
            fida_data=read_fidasim(filename,'eq','geom','spec','b3d',b3d_name);
        elseif  ~isempty(efit)
            fida_data=read_fidasim(filename,'eq','geom','spec','efit',efit);
        elseif  ~isempty(vmec)
            fida_data=read_fidasim(filename,'eq','geom','spec','vmec',vmec);
        end
    elseif lxerr
        fida_data=read_fidasim(filename,'geom','spec','eq');        
    else
        fida_data=read_fidasim(filename,'geom','spec');
    end
    if isfield(fida_data,'spec')
        data=fida_data.spec;
    end
    if isfield(fida_data,'geom')
        geom=fida_data.geom;
    end
%end
if isempty(data)
    disp([' Reading file: ' spec_name]);
    data=read_hdf5(spec_name);
end
if isempty(geom)
    disp([' Reading file: ' geom_name]);
    geom=read_hdf5(geom_name);
end

R = data.radius;
full =data.full;
half= data.half;
third= data.third;
halo= data.halo;
dcx= data.dcx;
fida= data.fida;
brems= data.brems;
lambda =data.lambda;

if isfield(data,'pfida')
    pfida=data.pfida;
    disp('Passive FIDA in Spectrum!')
else
    pfida=zeros(size(fida));
end

spec = full+ half + third + halo + dcx + fida + pfida;% + brems;



if lload_fidasim
    chan_description=channel;
    if ischar(channel)
        if strcmp(channel,'all')
            channel = true(size(deblank(geom.spec.id)));
        else
            channel = contains(deblank(geom.spec.id),channel);
        end
    elseif iscell(channel)
        channel_tmp = false(geom.spec.nchan,1);
        for i=1:numel(channel)
            channel_tmp(:,i) = contains(deblank(geom.spec.id),channel{i});
        end
        channel = channel_tmp;
    elseif isnumeric(channel)
        channel_tmp=false(geom.spec.nchan,1);
        if max(channel)>numel(channel_tmp)
            disp('Requested channel outside bounds of channels found in geometry file! Exiting!')
            return
        end
        channel_tmp(channel)=true;
        channel = channel_tmp;
    end
    dex=channel;
end

if isfield(in_data,'instfu_gamma')
    cwav_mid=mean(lambda);
    %cwav_mid = interp1(1:size(lambda,1),lambda,size(lambda,1)/2.);
    %if ~lload_fidasim
    %Flipud is necessary to emulate fplot.pro behavior from FIDASIM4
    instfu = flipud(box_gauss_funct(lambda,0.,1.,cwav_mid,in_data.instfu_gamma,in_data.instfu_box_nm));
    disp(['Applying Instrument function to FIDASIM data: ', filename]);
    if size(spec,2)~=size(instfu,2)
        disp('Careful! Only applying Instrument function to known LOS!')
    end
    for i = 1:size(instfu,2)
        spec(:,i) = conv(spec(:,i),instfu(:,i),'same');
    end
end
if lmean == 1
    spec = movmean(spec,15);
    disp(['Applying moving mean to FIDASIM data: ', filename]);
end

% dispersion_tmp = diff(lambda_dat,1);
dispersion_tmp = diff(lambda,1);
dispersion_tmp = [dispersion_tmp; dispersion_tmp(end,:)];
dispersion_tmp = repmat(dispersion_tmp,1,size(spec,2));

bg_dex = (lambda > bg_range(1)) & (lambda < bg_range(2));
bg = sum(brems.*dispersion_tmp.*bg_dex,1,'omitnan')./sum(dispersion_tmp.*bg_dex,1,'omitnan');

m=size(full',2);
[val,loc] = max(  fliplr(logical(full')),  [],2);
i=m+1-loc;
i(val==0)=m;
[val,loc] = max(  (logical(full')),  [],2);
k2=m+1-loc;
k2(val==0)=m;
bes_range_calc=[lambda(k2),lambda(i)];
if (size(bes_range,1)==numel(bg) && ~lload_fidasim) | lspecbes
    bes_dex = (lambda > repmat(bes_range(:,1)',size(lambda,1),1)) & (lambda < repmat(bes_range(:,2)',size(lambda,1),1));
    bes_dex=[bes_dex,zeros(size(bes_dex,1),numel(bg)-size(bes_dex,2))];
    % bes = sum(spec.*dispersion_tmp.*bes_dex,1,'omitnan');
    bes_range_tmp=bes_range;
    bes_range=zeros(size(bes_range_calc));
    bes_range(boolean(sum(bes_dex,1)),:)=bes_range_tmp(boolean(sum(bes_dex,1)),:);
    bes_range(~sum(bes_dex,1),:)=bes_range_calc(~sum(bes_dex,1),:);
    disp('BES from complete FIDASIM spectrum!')
    bes=trapz(lambda,spec.*bes_dex);
else
    disp('BES from full Energy component only!')
    bes = sum(full.*dispersion_tmp,1,'omitmissing')/3;%Approximate BES by full Beam component
    if isempty(bes_range)
        bes_range=bes_range_calc;
    else
        bes_range=[bes_range; repmat(fida_range,numel(bes)-size(bes_range,1),1)];
    end
end

fida_dex = (lambda > fida_range(1)) & (lambda < fida_range(2));
fida = sum(spec.*dispersion_tmp.*fida_dex,1,'omitnan');

if fac~=1
    dispname = ['', name, ', scaling factor: ' num2str(fac)];
else
    dispname = ['', name];
end


for i = 1:size(plot_type,2)
    if i>numel(ax)
        figure;
        ax{i} = gca;
        hold on;
    end
    if lrho==1 && isfield(geom.spec,'rho')
        [~,I]=sort(R(dex));
        R_plt=geom.spec.rho(dex);
        R_plt=R_plt(I);
        xlabel(ax{i},'\rho_{tor} [-]')
        xlim(ax{i},[0 1])
    else
        [R_plt,I]=sort(R(dex));
        xlabel(ax{i},'R [cm]')
    end
    switch lower(plot_type{i})
        case 'spectrum'
            plot_fidasim(filename,'spectrum',channel_spec,'ax',ax(i),'sim_data',in_data);
        case 'bck'
            tmp = bg(dex);
            ystr = 'BG [Ph/(s m^2 sr)]';
        case 'bes'
            tmp = bes(dex);
            ystr = 'BES [Ph/(s m^2 sr)]';
        case 'fida'
            tmp = fida(dex);
            ystr = 'FIDA [Ph/(s m^2 sr)]';
        case 'fidabes'
            tmp = fida(dex)./bes(dex).*abs(diff(bes_range(dex,:),1,2)'./diff(fida_range));
            ystr = 'FIDA/BES [-]';
        case {'fida_bck','fidabck','fidabg'}
            tmp = fida(dex)./bg(dex).*diff(bg_range)./diff(fida_range);
            ystr = 'FIDA/BACKGROUND [-]';
        case 'fidaspec'
            plot(ax{i},lambda(fida_dex), spec(fida_dex,dex),linestyle,'DisplayName',dispname, 'LineWidth',2.0);
            xlabel(ax{i},'Wavelength [nm]')
            ylabel(ax{i},'Intensity [Ph/(s nm m^2 sr)]')
            continue
        case 'bes_comp'
            %spec = full+ half + third + halo + dcx + fida + pfida;% + brems;
            bestmp(:,1)=  sum(full.*dispersion_tmp.*bes_dex,1,'omitnan');
            bestmp(:,2)=  sum(half.*dispersion_tmp.*bes_dex,1,'omitnan');
            bestmp(:,3)=  sum(third.*dispersion_tmp.*bes_dex,1,'omitnan');
            bestmp(:,4)=  sum((halo+dcx).*dispersion_tmp.*bes_dex,1,'omitnan');
            bestmp(:,5)=  sum(fida.*dispersion_tmp.*bes_dex,1,'omitnan');
            bestmp(:,6)=  sum(pfida.*dispersion_tmp.*bes_dex,1,'omitnan');
            bar(ax{i},R_plt,bestmp(dex,:),'stacked');
            continue


    end

    if fac~=1
        tmp = tmp.*fac;
    end
    if liota
        [pts(1,:),pts(2,:),pts(3,:)]=getPointsAlongAxis(geom.nbi.src+geom.nbi.axis*500, geom.nbi.axis,1800,300,true);
        if ~isempty(efit)
            psi_line=interp2(efit.xgrid*100,efit.zgrid*100,efit.psixz',pts(1,:),pts(3,:),'linear',NaN);
            q_line=interp1(linspace(efit.psiaxis,efit.psilim,numel(efit.qpsi)),efit.qpsi,psi_line);
            for m=1:numel(iota_vec)
                [~,tmpp]=findpeaks(-abs(q_line-iota_vec(m)));
                if ~isempty(tmpp)
                    x_plt{m}=tmpp;
                else
                    disp(['Found no surfaces for q=', num2str(iota_vec(m))])
                end
            end
            xline(ax{i},efit.xax(1)*100*fac,'--',{sprintf('Axis',iota_vec(m))},'HandleVisibility','off');
        elseif ~isempty(vmec)
            if ~isfield(vmec,'Fchi')
                vmec=vmec_rzphi_s_interp(vmec);
            end
            phin = vmec.phi./vmec.phi(end);
            s_line = vmec.Fchi(pts(1,:)/100,mod(pts(2,:),vmec.zeta(end)),pts(3,:)/100);
            q_line=interp1(phin,vmec.iotaf,s_line);
            x_plt={};
            for m=1:numel(iota_vec)
                [~,tmpp]=findpeaks(-abs(q_line-iota_vec(m)),'MinPeakHeight',-.01);
                if ~isempty(tmpp)
                    x_plt{m}=tmpp;
                else
                    disp(['Found no surfaces for q=', num2str(iota_vec(m))])
                end
            end
            %xl=xline(ax{i},vmec.axisr,'--',{sprintf('Axis',iota_vec(m))},'HandleVisibility','off');
        end
        if lrho
            if ~isempty(vmec)
                pts(1,:) = sqrt(vmec.Fchi(pts(1,:)/100,mod(pts(2,:),vmec.zeta(end)),pts(3,:)/100));
            else
                pts(1,:) = interp3(fida_data.eq.fields.r,fida_data.eq.fields.phi,fida_data.eq.fields.z,...
                    permute(sqrt(fida_data.eq.fields.s),[ 3 1 2]),pts(1,:),mod(pts(2,:),fida_data.eq.fields.phi(end)),pts(3,:),'linear',NaN);
            end
        end
        for m=1:numel(iota_vec)
            xl=xline(ax{i},pts(1,x_plt{m}),'-',{sprintf('q=%.2f',iota_vec(m))},'HandleVisibility','off');
            for j=1:numel(xl)
                xl(j).LabelHorizontalAlignment='center';
                xl(j).LabelVerticalAlignment='bottom';
            end
        end
    end

    tmp=tmp(I);
    R_plt=reshape(R_plt,[1 numel(R_plt)]);
    if ~lxerr
        plot(ax{i},R_plt, tmp,linestyle,'DisplayName',dispname, 'LineWidth',2.0);
        %plot(ax{i},R_plt(channel), tmp(channel_spec),'.','DisplayName',['Channel ',num2str(channel)], 'LineWidth',2.0);
    else
        neut=read_fidasim(filename,'neut').neut;
        input=fida_data.input;
        eq=fida_data.eq;
        [upts,vpts,wpts]=getPointsAlongAxis(geom.spec.lens, geom.spec.axis,10000,400,false);
        uvwpts=[upts(:),vpts(:),wpts(:)];
        [phipts,rpts,zpts]=cart2pol(upts/100,vpts/100,wpts/100);
        %denf_pts=interp3(eq.fields.r,eq.fields.z,eq.fields.phi,permute(denf,[2 1 3]),rpts*100,zpts*100,mod(phipts,eq.fields.phi(end)));
        neut.dens = neut.fdens+neut.hdens+neut.tdens+neut.dcxdens+neut.halodens;
        neut.dens=squeeze(sum(neut.dens,1));%Sum over all levels
        ndisc=256;
        if lrho
            if numel(b3d_name)~=0
                rho_pts=real(sqrt(beams3d_getvals(b3d_name,rpts,phipts,zpts,'S_ARR').S_ARR));
            elseif  ~isempty(efit)
                rho_pts=interp2(efit.xgrid,efit.zgrid,(efit.psixz'-efit.psiaxis)./(efit.psilim-efit.psiaxis),rpts,zpts,'linear',NaN);
            elseif  ~isempty(vmec)
                rho_pts = vmec.Fchi(rpts,mod(pts(2,:),vmec.zeta(end)),zpts);
            end
            xedges=linspace(0,1.5,ndisc);
            xplt=xedges+(xedges(2)-xedges(1))/2;
            [discr,~]=discretize(rho_pts,xedges);
            discr(isnan(discr))=ndisc+1;
        else
            %xedges=linspace(min(rpts(~isnan(dens_pts)),[],'all','omitnan'),max(rpts(~isnan(dens_pts)),[],'all','omitnan'),ndisc);
            xedges=linspace(eq.fields.r(1),eq.fields.r(end),ndisc);
            xplt=(xedges+(xedges(2)-xedges(1))/2);
            [discr,~]=discretize(rpts*100,xedges);
            discr(isnan(discr))=ndisc+1;
        end

        %Densi in beam grid coordinates
        %UVW Corrdinates are Machine coordinates, XYZ are beam grid coordinates
        %coords= xyz_to_uvw(input.alpha, input.beta, input.gamma, coords, input.origin);
        [xyz] = uvw_to_xyz(input.alpha, input.beta, input.gamma, uvwpts, input.origin);
        xyzpts=reshape(xyz,size(rpts,1),size(rpts,2),3);
        xpts=xyzpts(:,:,1);
        ypts=xyzpts(:,:,2);
        zpts=xyzpts(:,:,3);
        [xg,yg,zg]=ndgrid(neut.grid.x,neut.grid.y,neut.grid.z);
        dens=griddedInterpolant(xg,yg,zg,neut.dens,'linear','none');
        dens_pts=dens(xpts,ypts,zpts);
        los_ind=repmat(1:geom.spec.nchan,size(xpts,1),1);

        xtmp=dens_pts(:);%.*denf_pts(:);
        xtmp = accumarray([discr(:),los_ind(:)],xtmp,[ndisc+1 geom.spec.nchan+1],@(x) sum(x,'omitmissing'));
        xtmp=xtmp(1:end-1,1:end-1);
        %Normalize
        xtmp=xtmp./sum(xtmp,1,'omitmissing'); %Normalize each channel

        %Half height
        hheight=(max(xtmp,[],1,'omitmissing')-min(xtmp,[],1,'omitmissing'))/2;
        index1=ones(1,size(xtmp,2));
        index2=ones(1,size(xtmp,2));
        for j=1:size(xtmp,2)
            if ~isnan(hheight(j))
            index1(j)=find(xtmp(:,j)>=hheight(j),1,'first');
            index2(j)=find(xtmp(:,j)>=hheight(j),1,'last');
            end
        end
        neg=xplt(index1(dex));
        pos=xplt(index2(dex));
        
        neg=R_plt-neg(I);
        pos=pos(I)-R_plt;
        errorbar(ax{i},R_plt, tmp,neg,pos,'horizontal',linestyle,'DisplayName',dispname, 'LineWidth',2.0);

    end
    ylabel(ax{i},ystr)

    if lsave
        sname = [name, '_', plot_type{i}];
        if lrho
            sname = [sname,'_rho'];
        end
        savefig(gcf,[sname,'.fig']);
        if leps
            exportgraphics(ax{i}.Parent,[sname,'.eps'],'Resolution',600,'BackgroundColor','none');
        end
        exportgraphics(ax{i}.Parent,[sname,'.png'],'Resolution',600);
    end
end
plot_data.R = R;
plot_data.R_plt = R_plt;
plot_data.R_plt_ind=find(dex);
plot_data.bg = bg;
plot_data.bes = bes;
plot_data.fida = fida;
plot_data.fidabes = fida./bes;
plot_data.spec = spec;
plot_data.lambda=lambda;
plot_data.bg_range=bg_range;
plot_data.bes_range=bes_range;
plot_data.fida_range=fida_range;
plot_data.ax=ax;
plot_data.dex = dex;
end

% function F = box_gauss_funct(X,A,B,C,D,E) % From /afs/ipp/home/s/sprd/XXX_DIAG/LIB
% gam   = double(D);
% width = double(E);
% rl    = abs(0.5d0*width./gam);
% Z     = abs((double(X)-double(C))./gam);
% F     = double(B)*(0.5d0./width.*(erf(Z+rl) - erf(Z-rl)))+double(A);
%
% % Normalization and cutoff
% %F = F./sum(F,1);
% %F(F<1e-5) = 0;
% end