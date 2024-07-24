function plot_data = plot_fidasim_profiles(filename,in_data,varargin)
%PLOT_FIDASIM_PROFILES plots radial profiles from the FIDASIM code. As
%these are used to compare to experimental data, the integration range and
%dispersion is set by variables in the in_data structure generated e.g. by
%the get_bes_fida_aug_data function. Plots are shared between the
%two functions.
%
% Example usage
%      [~,in_data] = get_bes_fida_aug_data(filename,'t_point',3.5,'fidabes');
%      plot_fidasim_profiles(filename,in_data,'X);
%      !!! 'X' can be 'fida', 'bes', or 'fidabes'
%
% Miscellaneous Arguments
%      plot_fidasim(runid,'mean'); %Apply moving mean to spectrum
%      plot_fidasim(runid,'in_data',in_data); %Used for dispersion
%      plot_fidasim(runid,'save'); %Export figures (.fig and .png)
%      plot_fidasim(runid,'name', 'test'); %ID Name for legend
%      plot_fidasim(runid,'fac', 1.0); %Scaling factor
%
bg_range=in_data.bg_range;

bes_range = in_data.bes_range;
fida_range = in_data.fida_range;
dispersion = in_data.dispersion;
lambda_dat = in_data.lambda;
if isfield(in_data,'dex')
    dex = in_data.dex;
else
    dex = [];
end


if isfield(in_data,'ax')
    ax=in_data.ax;
else
    ax = {};
end

lsave = 0;
lmean = 0;
lload_fidasim=0;
lrho=0;
leps=0;
plot_type = {};
linestyle = '+';
fac = 1;

if strcmp(filename(end-10:end),'_spectra.h5')
    filename=filename(1:end-11);
end
name = filename(1:end-2);

if nargin > 2
    i = 1;
    while i < nargin-1
        switch varargin{i}
            case {'FIDA','BES','FIDABES','fida','bes','fidabes','bck'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
            case 'mean'
                lmean =1;
            case 'save'
                lsave = 1;
            case 'avg_frames'
                i = i+1;
                avg_frames = varargin{i};
            case 'rho'
                lrho=1;
                i=i+1;
                b3d_name=varargin{i};
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



spec_name = [filename,'_spectra.h5'];
geom_name = [filename,'_geometry.h5'];

if ~isfile(geom_name)|| lload_fidasim
    if lrho
        fida_data=read_fidasim(filename,'eq','geom','spec','b3d',b3d_name);
    else
        fida_data=read_fidasim(filename,'geom','spec');
    end
    data=fida_data.spec;
    geom=fida_data.geom;
else
    data=read_hdf5(spec_name);
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
else
    pfida=zeros(size(fida));
end

spec = full + half + third + halo + dcx + fida + pfida;% + brems;



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

    end
    dex=channel;
end


cwav_mid=mean(lambda);
if ~lload_fidasim
instfu = box_gauss_funct(lambda,0.,1.,cwav_mid,in_data.instfu_gamma,in_data.instfu_box_nm);
disp(['Applying Instrument function to FIDASIM data: ', filename]);
if numel(in_data.instfu_gamma)==numel(in_data.names)
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

if size(bes_range,1)==numel(bg)&&~lload_fidasim
bes_dex = (lambda > repmat(bes_range(:,1)',size(lambda,1),1)) & (lambda < repmat(bes_range(:,2)',size(lambda,1),1));
bes = sum(spec.*dispersion_tmp.*bes_dex,1,'omitnan');
else
    bes = sum(full.*dispersion_tmp,1,'omitnan')/3;%Approximate BES by full Beam component
end

fida_dex = (lambda > fida_range(1)) & (lambda < fida_range(2));
fida = sum(spec.*dispersion_tmp.*fida_dex,1,'omitnan');

for i = 1:size(plot_type,2)
    if i>numel(ax)
        figure;
        ax{i} = gca;
        hold on;
    end
    switch lower(plot_type{i})
        case 'bck'
            tmp = bg(dex);
            ystr = 'BG';
        case 'bes'
            tmp = bes(dex);
            ystr = 'BES';
        case 'fida'
            tmp = fida(dex);
            ystr = 'FIDA';
        case 'fidabes'
            tmp = fida(dex)./bes(dex);
            ystr = 'FIDA/BES';
            % if lsave
            %     legend(ax{i},'Location','southwest');
            % end

    end

    if fac~=1
        tmp = tmp.*fac;
        dispname = ['', name, ', scaling factor: ' num2str(fac)];
    else
        dispname = ['', name];
    end
    if lrho==1 && isfield(geom.spec,'rho')
        plot(ax{i},geom.spec.rho(dex), tmp,linestyle,'DisplayName',dispname, 'LineWidth',2.0);
        xlabel(ax{i},'\rho_{tor} [-]')   
        xlim([0 1])
    else
        plot(ax{i},R(dex), tmp,linestyle,'DisplayName',dispname, 'LineWidth',2.0);
        xlabel(ax{i},'R [cm]')
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
plot_data.bes = bes;
plot_data.fida = fida;
plot_data.spec = spec;
plot_data.lambda=lambda;
plot_data.ax=ax;
plot_data.dex = dex;

end


function F = box_gauss_funct(X,A,B,C,D,E) % From /afs/ipp/home/s/sprd/XXX_DIAG/LIB
gam   = double(D);
width = double(E);
rl    = abs(0.5d0*width./gam);
Z     = abs((double(X)-double(C))./gam);
F     = double(B)*(0.5d0./width.*(erf((Z+rl)) - erf((Z-rl))))+double(A);

% Normalization and cutoff
F = F./sum(F,1);
F(F<1e-5) = 0;
end