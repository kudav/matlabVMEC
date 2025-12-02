function [ ax, plt_data ] = plot_fidasim(file,varargin)
%PLOT_FIDASIM Makes plots of FIDASIM data generated with the BEAMS3D
%Interface.
%The PLOT_FIDASIM function creates various canned plots of the files used
%for running FIDASIM and the outputs. The function reads the necessary
%files automatically when supplied a FIDASIM runid. The function is quite
%flexible and can also be used to compare data from multiple runs and
%across codes. For plotting FIDA/BES profiles, see plot_fidasim_profiles.
%
% Example usage
%      plot_fidasim(f); %f can be id or struct from read_fidasim
%      plot_fidasim(f,'overview'); %NOT IMPLEMENTED
%      plot_fidasim(f,'profiles'); %Kinetic profiles
%      plot_fidasim(f,'spectrum', channel_no); %Spectral components
%      plot_fidasim(f,'fslice'); %Outboard midplane profile of dist.(f)
%      plot_fidasim(f,'ba'); %All magnetic field components on midplane
%      plot_fidasim(f,'bx_'); %Magnetic field (x=r,p,z)
%      plot_fidasim(f,'denf_'); %FI density on RZ plane from denf
%      plot_fidasim(f,'fdenf_'); %FI density on RZ plane from f
%       optionally followed by index for toroidal/vertical position
%       '_' can be:     none: profile@midplane, index for tor. pos
%                       '2d': R-Z plane
%                           'torint' calculates toroidal integral
%                           'intersection' calculates LOS positions
%                           'sep' shows approx. separatrix pos.
%                       'tor': midplane (R-PHI)
%
%      plot_fidasim(f,'energy'); %Energy dist, int. over real space
%      plot_fidasim(f,'pitch'); %Pitch dist, int. over real space
%      plot_fidasim(f,'epplot'); %E-p dist.,  int. over real space
%      plot_fidasim(f,'epplot',[R PHI Z]); %E-p dist., at spatial point
%      plot_fidasim(f,'epplot',[R1 R2 P1 P2 Z1 Z2]); %E-p int over
%                               spatial range
%      plot_fidasim(f,'bmir'); %Mirror field dist. (broken)
%      plot_fidasim(f,'q2d'); %Approx. of safety factor
%      plot_fidasim(f,'ndensvert'); %Neutral density, vertical
%      plot_fidasim(f,'ndenshorz'); %Neutral density, horizontal
%      plot_fidasim(f,'ndenscross'); %Neutral density, cross (beam normal)
%                      replace 'ndens' with 'halo' for dcx+halo density
%      plot_fidasim(f,'weights',[lamda,channel_num]); %FIDA weights
%      plot_fidasim(f,'weights'); %FIDA weights
%      plot_fidasim(f,'ax', ax); %Figure handles for sharing plots
%
% Miscellaneous Arguments
%      plot_fidasim(f,'mean'); %Apply moving mean to spectrum
%      plot_fidasim(f,'contour', levels); %Show area as contour plot
%      plot_fidasim(f,'style', linestyle); %Define linestyle
%      plot_fidasim(f,'sim_data',sim_data); %Used for dispersion
%      plot_fidasim(f,'save'); %Export figures (.fig and .png)
%      plot_fidasim(f,'name', 'test'); %ID Name for legend
%      plot_fidasim(f,'fac', 1.0); %Scaling factor
%
% Example usage of comparing BEAMS3D and TRANSP results:
%     filename_b3d = fidasim_b3d
%     filename_transp = transp
%     [ax, n_b3d] = plot_fidasim(filename_b3d,'energy','pitch', 'profiles');
%     [ax, n_transp] = plot_fidasim(filename_transp,'energy','pitch','profiles','ax',ax);
%
% Maintained by: David Kulla (david.kulla@ipp.mpg.de)
% Version:       1.00

ec=1.6021773300E-19; % Charge of an electron (leave alone)
amu = 1.66053906660E-27; % Dalton [kg]


n_fida=-1;
input={};
eq={};
dist={};
geom={};
rreq_ind=0;
zreq_ind=0;
phireq_ind=0;

if ischar(file)
    if (strcmp(file(end-2:end),'.h5'))
        disp('ERROR: Only give runid (dist and eq files are loaded automatically!');
        disp(['       Filename: ' file]);
    end
    dist_name = [file, '_distribution.h5'];
    eq_name = [file,'_equilibrium.h5'];
    neut_name = [file, '_neutrals.h5'];
    spec_name = [file,'_spectra.h5'];
    geom_name = [file,'_geometry.h5'];
    birth_name = [file,'_birth.h5'];
    weight_name = [file,'_fida_weights.h5'];
    lloaded=0;
    name = file;
else
    lloaded=1;
    if isfield(file,'eq')
        eq=file.eq;
    end
    if isfield(file,'dist')
        dist=file.dist;
        dr = dist.r(2)-dist.r(1);
        dz = dist.z(2)-dist.z(1);
        nr = double(dist.nr);
        nz = double(dist.nz);
        if ndims(dist.f) == 5
            dphi = dist.phi(2) - dist.phi(1);
            rdphi= dist.r*dphi;
            nphi=double(dist.nphi);
            area=dr.*dz;
            % Volume (function of R)
            vol = dist.r.*dphi.*area;
            vol2d=repmat(vol,[1 dist.nz dist.nphi]);
            n_fida = sum(dist.denf.*vol2d,'all');
        else
            dphi=2*pi;
            rdphi=dist.r*dphi;
            nphi=1;
            n_fida = 2*pi*dr*dz*sum(dist.r.*sum(squeeze(dist.denf(:,:,1)),2)); %Axisymmetric only.
        end

        [~,zreq_ind]=min(abs(dist.z));
    end
    if isfield(file,'geom')
        geom=file.geom;
    end
    if isfield(file,'neut')
        neut=file.neut;
    end
    if isfield(file,'weight')
        weight=file.weight;
    end
    if isfield(file,'spec')
        spec=file.spec;
    end
    if isfield(file,'birth')
        birth=file.birth;
    end
    if isfield(file,'input')
        input=file.input;
        name=input.runid;
        file=input.runid;
    else
        name='FIDASIM';
        file='LOADED FIDASIM';
    end
end



lsave = 0;
plot_type = {};
ax = {};
fac = 1;

leq = 0;
ldist = 0;
lspec = 0;
lneut = 0;%
lmean=0;% Mean over volume instead of integration
lgeom =0; %Load geometry
lweight=0; %Load weight
lbirth=0; % Load birth
lcontour=0; %Plot contour
levels=2; %Number of levels for contour
linput=1; % read input namelist
ltorint=0; %Integrate over toroidal direction
ltor=0; %Plot Toroidal cut
lz=0; %Whether to plot horizontal cut
zval=0.0; %Requested z value
liota=0;%to plot iota values
efit={};
vmec={};
lintersection=0; %Plot of cutplane intersection with LOS
sim_data = {}; %Machine-Specific data
channel = 0; %Requested channels
linestyle = '-';
index=1;
index_in=[];
rotation=0;
lpassive=0;
lbrems=0;
llegend=0;
leps=0;
lsep=0;
ltrim=0; %Trim first and last toroidal gridpoints for distribution
ldiff=0; %Plot difference to second distribution
lrel=0;
dist2={};
dist2_name='';
tmp=[];
tmp2=[];
plot_type={};
if nargin > 1
    i = 1;
    while i < nargin
        switch lower(varargin{i})
            case {'overview','profiles','profiles_rho',...
                    'ba','br2d','bt2d','bz2d',...
                    'brtor','bttor','bztor','q2d',...
                    'te2d','ne2d','ti2d', 'er2d','et2d','ez2d',...
                    'vt2d','zeff2d','denn2d','mask2d'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                leq = 1;
                if numel(varargin)>i
                    if ~ischar(varargin{i+1})
                        i=i+1;
                        index_in = varargin{i};
                    end
                end

            case {'ep2d'}
                varargin{i}='epplot';
                continue;
            case {'fslice','denf2d','denf','denftor',...
                    'fdenf','fdenf2d','fdenftor','fdenf3d',...
                    'pitch','energy','epplot',...
                    'bmir'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                ldist = 1;
                leq=1;
                if numel(varargin)>i
                    if ~ischar(varargin{i+1})
                        i=i+1;
                        index_in = varargin{i};
                    end
                end
            case 'vflow2d'
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                ldist=1;
                leq=1;
                linput=1;
                if numel(varargin)>i
                    if ~ischar(varargin{i+1})
                        i=i+1;
                        index_in = varargin{i};
                    end
                end
            case {'lcfs','sep','separatrix'}
                lsep=1;
                leq=1;
            case {'trim'}
                ltrim=1;
            case {'weights','weight_dist','weights_dist'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                if strcmp(varargin{i}(end-2:end),'ist')
                    ldist = 1;
                end
                lweight=1;
                lgeom=1;
                if numel(varargin)>i
                    if ~ischar(varargin{i+1})
                        i=i+1;
                        index_in = varargin{i};
                    end
                end
            case{'ndensvert', 'ndenshorz', 'ndenscross'...
                    'fdensvert', 'fdenshorz', 'fdenscross'...
                    'halovert', 'halohorz', 'halocross'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                lneut = 1;
                lgeom=1;
                linput=1;
                if numel(varargin)>i
                    if ~ischar(varargin{i+1})
                        i=i+1;
                        index_in = varargin{i};
                    else
                        index_in =20;
                    end
                end
            case{'ndens2d', 'ndenstor','ndens'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                lneut = 1;
                lgeom=1;
                linput=1;
                leq=1;
                if numel(varargin)>i
                    if ~ischar(varargin{i+1})
                        i=i+1;
                        index_in = varargin{i};
                    else
                        index_in = 1;
                    end
                end
                index=index_in;
            case{'spectrum'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                lspec = 1;
                lgeom = 1;
                linput=1;
                i=i+1;
                channel = varargin{i};
            case{'los3d', 'lostor','los2d'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                %lspec = 1;
                lgeom = 1;
                linput=1;
                i=i+1;
                channel = varargin{i};
                i=i+1;
                length=varargin{i};
                %                 i=i+1;
                %                 color=varargin{i};
            case {'fida', 'bes', 'fidabes'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                disp(['ERROR: Option ', varargin{i}, ' not implemented here. Use plot_fidasim_profiles instead.']);
                lspec = 1;
                lgeom = 1;
            case {'birth_r','birth_z','birth_pitch', 'birth_phi'...
                    'birth_r_gc','birth_z_gc','birth_phi_gc' }
                lbirth = 1;
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
            case {'diff','diffrel','diff_rel'}
                ldiff=1;
                 if strcmp(varargin{i}(end-2:end),'rel')
                     lrel=1;
                 end
                i=i+1;
                file2=varargin{i};
                dist2_name = [file2, '_distribution.h5'];
            case 'channel'
                i=i+1;
                channel=varargin{i};
            case 'eps'
                leps = 1;
            case 'legend'
                llegend = 1;
            case 'mean'
                lmean = 1;
            case {'contour','contours'}
                lcontour=1;
                if numel(varargin)>i
                    if ~ischar(varargin{i+1})
                        i=i+1;
                        levels = varargin{i};
                    end
                end
            case 'sim_data'
                i=i+1;
                sim_data = varargin{i};
            case 'intersection'
                lintersection=1;
                lgeom=1;
            case 'qeqdsk'
                liota=1;
                lgeom=1;
                i=i+1;
                efit=read_efit(varargin{i});
                i=i+1;
                iota_vec=varargin{i};
            case 'iotavmec'
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
            case 'save'
                lsave = 1;
            case 'passive'
                lpassive=1;
            case 'brems'
                lbrems=1;
            case 'torint'
                ltorint=1;
            case 'ax'
                i=i+1;
                ax = varargin{i};
            case 'fac'
                i = i+1;
                fac = varargin{i};
            case 'name'
                i = i+1;
                name = varargin{i};
            case 'style'
                i = i+1;
                linestyle = varargin{i};
            case 'noinput'
                linput=0;
            case {'frominputs','frominput','fromdat'}
                linput=1;                
            case 'z'
                i=i+1;
                zval=varargin{i};
            otherwise
                disp(['ERROR: Option ', varargin{i}, ' not found!']);
        end
        if numel(plot_type{end})>2
            if strcmp(plot_type{end}(end-1:end),'2d')
                ltor=1;
            elseif strcmp(plot_type{end}(end-2:end),'tor')
                lz=1;
            end
        end
        i=i+1;
    end
end

if ltor&&lz
    disp('Request either toroidal or horizontal cut, not both!')
    return
end
% if ~isempty(index_in)
%    index = index_in(end);
% end
if ~lloaded
    if linput
        if isfile([file,'_inputs.dat'])
        input=read_namelist([file,'_inputs.dat'],'fidasim_inputs');
        else
            disp('Input namelist not found! Using standard names!')
            input={};
            linput=false;
        end
    end
    if ldist
        if linput && isfile(input.distribution_file)
            dist_name=input.distribution_file;
        elseif isfile(dist_name)
            dist_name=dist_name;
        else
            disp('ERROR: Distribution file not found, check filename!');
            disp(['       Filename: ' dist_name]);
            disp(['       Filename: ' input.distribution_file]);
            return
        end
        groupnames={'energy','pitch','r','z','phi','nenergy','npitch','nr','nz','nphi'};
        for k=1:numel(groupnames)
            dist.(groupnames{k})= h5read(dist_name,['/',groupnames{k}]);
        end
        %dist = read_hdf5(dist_name);
    end
    if leq
        if isfile(eq_name)
            eq = read_hdf5(eq_name);
        elseif linput && isfile(input.equilibrium_file)
            eq = read_hdf5(input.equilibrium_file);
        else
            disp('ERROR: Equilbirium file not found, check filename!');
            disp(['       Filename: ' input.equilibrium_file]);
            if ldist
                eq={};
                eq.fields.z=dist.z;
                eq.fields.r=dist.r;
            else
                disp('No R/Z information available, exiting!')
                return
            end
        end
        [~,zreq_ind]=min(abs(eq.fields.z-zval));
    end
end
if leq
    %Assign dimension slice indices
    [phi0_ind,r0_ind,z0_ind,e_min,p_min]=deal(1);
    phi1_ind=eq.fields.nphi;
    r1_ind=eq.fields.nr;
    z1_ind=eq.fields.nz;
end
if ldist
    e_max=dist.nenergy;
    p_max=dist.npitch;
end
if isscalar(index_in)
    if ltor
        phireq_ind=index_in;%
    elseif lz
        zreq_ind=index_in;%
    end
elseif numel(index_in)==3
    [~,rreq_ind]=min(abs(eq.fields.r-index_in(1)));
    [~,phireq_ind]=min(abs(eq.fields.phi-index_in(2)));
    [~,zreq_ind]=min(abs(eq.fields.z-index_in(3)));
elseif numel(index_in)==4
    [~,e_min]=min(abs(dist.energy-index_in(1)));
    [~,e_max]=min(abs(dist.energy-index_in(2)));
    [~,p_min]=min(abs(dist.pitch-index_in(3)));
    [~,p_max]=min(abs(dist.pitch-index_in(4)));
    dist.energy=dist.energy(e_min:e_max);
    dist.pitch=dist.pitch(p_min:p_max);
elseif numel(index_in)==5
    [~,e_min]=min(abs(dist.energy-index_in(1)));
    [~,e_max]=min(abs(dist.energy-index_in(2)));
    [~,p_min]=min(abs(dist.pitch-index_in(3)));
    [~,p_max]=min(abs(dist.pitch-index_in(4)));
    if ltor
        phireq_ind=index_in(5);%
    elseif lz
        zreq_ind=index_in(5);
    end
elseif numel(index_in)==6
    [~,r0_ind]=min(abs(eq.fields.r-index_in(1)));
    [~,z0_ind]=min(abs(eq.fields.z-index_in(5)));
    [~,r1_ind]=min(abs(eq.fields.r-index_in(2)));
    [~,z1_ind]=min(abs(eq.fields.z-index_in(6)));
    if eq.fields.nphi > 1
        [~,phi0_ind]=min(abs(eq.fields.phi-index_in(3)));
        [~,phi1_ind]=min(abs(eq.fields.phi-index_in(4)));
    end

end


if ~lloaded
    if ldist
        if ldiff
            groupnames={'energy','pitch','r','z','phi','nenergy','npitch','nr','nz','nphi'};
            for k=1:numel(groupnames)
                dist2.(groupnames{k})= h5read(dist2_name,['/',groupnames{k}]);
            end
            [dist2, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = load_dist_and_ind(index_in, dist2, dist2_name, phireq_ind, eq, rreq_ind, z0_ind, zreq_ind, e_min, p_min, e_max, p_max, ltor, phi1_ind, lz, zval, ltrim, r0_ind, r1_ind, z1_ind, phi0_ind);
        end
        [dist, phireq_ind, eq, z0_ind, zreq_ind, phi1_ind, zval, r0_ind, r1_ind, z1_ind, phi0_ind] = load_dist_and_ind(index_in, dist, dist_name, phireq_ind, eq, rreq_ind, z0_ind, zreq_ind, e_min, p_min, e_max, p_max, ltor, phi1_ind, lz, zval, ltrim, r0_ind, r1_ind, z1_ind, phi0_ind);
        % if ldiff
        %     dist.f=dist.f-dist2.f;
        %     dist.denf=dist.denf-dist2.denf;
        %     if lrel
        %         eps=1;%12e8;
        %         dist2.f(dist2.f<=eps)=eps;
        %         dist2.denf(dist2.denf==eps)=eps;
        %         dist2.f(isnan(dist2.f))=eps;
        %         dist2.denf(isnan(dist2.denf))=eps;  
        % 
        %         % dist.f=dist.f./dist2.f;                 
        %         % dist.denf=dist.denf./dist2.denf;
        % 
        %         dist.f(dist.f==eps)=eps;
        %         dist.denf(dist2.denf==eps)=eps;
        %         dist.f(isnan(dist.f))=eps;
        %         dist.denf(isnan(dist.denf))=eps;
        %     end
        % end
        dr = dist.r(2)-dist.r(1);
        dz = dist.z(2)-dist.z(1);
        nr = double(dist.nr);
        nz = double(dist.nz);    
        if dist.nphi>1||eq.fields.nphi>1
            dphi = dist.phi(2) - dist.phi(1);
            rdphi=dist.r*dphi;
            nphi=double(dist.nphi);
            area=dr.*dz;
            % Volume (function of R)
            vol = dist.r.*dphi.*area;
            vol2d=repmat(vol,[1 dist.nz dist.nphi]);
            n_fida = sum(dist.denf.*vol2d(r0_ind:r1_ind,z0_ind:z1_ind,phi0_ind:phi1_ind),'all');
        else
            dphi=2*pi;
            rdphi=dist.r*dphi;
            nphi=1;
            n_fida = 2*pi*dr*dz*sum(dist.r(r0_ind:r1_ind).*sum(squeeze(dist.denf(:,:,1)),2)); %Axisymmetric only.
        end
        if ltorint
            if nphi>1
                dist.f = trapz(dphi*nphi/(nphi-1),dist.r2d.*dist.f,5);
                dist.denf = trapz(dphi*nphi/(nphi-1),dist.r2d.*dist.denf,3);                
            else
                %dist.f=dist.r2d.*dist.f*2*pi;
                
            end
            index=1;
            cstring(end-2)='2';%Denote area density
        end    
        if zval~=dist.z(zreq_ind)
            [~,zreq_ind]=min(abs(dist.z-zval));
        end
    end

    if lweight
        weight= read_hdf5(weight_name);
    end
    if lneut
        neut = read_hdf5(neut_name);
        if ~isstruct(neut)
            disp('ERROR: Neutrals file not found, check filename!');
            disp(['       Filename: ' file]);
        end
    end

    if lspec
        spec = read_hdf5(spec_name);
        if ~isstruct(spec)
            disp('ERROR: Spectra file not found, check filename!');
            disp(['       Filename: ' file]);
        end
        %[~,I] = sort(spec.radius);
    end
    if lgeom
        if isfile(geom_name)
            geom = read_hdf5(geom_name);
        elseif linput && isfile(input.geometry_file)
            geom = read_hdf5(input.geometry_file);
        else
            disp('ERROR: Geometry file not found, check filename!');
            disp(['       Filename: ' geom_name]);
            disp(['       Filename: ' input.geometry_file]);
            lgeom=0;
        end
    end
    if lbirth
        birth = read_hdf5(birth_name);
        if ~isstruct(birth)
            disp('ERROR: Birth file not found, check filename!');
            disp(['       Filename: ' file]);
            lgeom=0;
        end
    end
else
    if leq
        [~,zreq_ind]=min(abs(eq.fields.z-zval));
    elseif ldist
        [~,zreq_ind]=min(abs(dist.z-zval));
    end
end

if lneut
    neut.dens = neut.fdens+neut.hdens+neut.tdens;%+neut.dcxdens+neut.halodens;
    neut.dens=squeeze(sum(neut.dens,1));%Sum over all levels
    neut.grid.vol = repmat(mean(diff(neut.grid.x))*mean(diff(neut.grid.y))*mean(diff(neut.grid.z)),size(neut.grid.x_grid));
    neut.nparts=neut.dens.*neut.grid.vol;
    nneutrals=1.d6*input.pinj/ (1.d3*input.einj*ec...
        *( input.current_fractions(1)      ...
        +  input.current_fractions(2)/2.d0 ...
        +  input.current_fractions(3)/3.d0 ) );
    if index_in==1
        [~,index]=min(abs(neut.grid.z));
    else
        index=index_in;
    end
end

if ischar(channel)
    chan_description = channel;
    if ~isempty(sim_data)
        channel = find(strcmp(channel,cellstr(deblank(sim_data.names'))));
        %channel = I(channel);%-1;
    elseif strcmp(channel,'all')
        channel = true(size(deblank(geom.spec.id)));
    else
        channel = contains(deblank(geom.spec.id),channel);
    end
elseif iscell(channel)
    chan_description=channel;
    if ~isempty(sim_data)
        channel_tmp = [];
        for i=1:numel(channel)
            channel_tmp =[channel_tmp, find(strcmp(channel{i},cellstr(deblank(sim_data.names'))))];
        end
        channel=sort(channel_tmp);
        %channel = I(channel);%-1;
    else
        channel_tmp = false(geom.spec.nchan,1);
        for i=1:numel(channel)
            %channel_tmp(:,i) = or(channel_tmp,contains(deblank(geom.spec.id),channel{i}));
            channel_tmp(:,i) = contains(deblank(geom.spec.id),channel{i});
        end
        channel = channel_tmp;
    end
elseif lgeom && channel~=0
    chan_description=geom.spec.id{channel};
    channel_tmp=channel;
    channel=false(geom.spec.nchan,1);
    channel(channel_tmp)=true;
end

if index==1
    index=index_in;
end

for i = 1:size(plot_type,2)
    if i>numel(ax)
        %figs{i}=figure;
        figure;
        ax{i} = gca;
        hold(ax{i},'on');
    else
        %allAxesInFigure = findall(figs{i},'type','axes');
        %ax{i} = allAxesInFigure(~ismember(get(allAxesInFigure,'Tag'),{'legend','Colobar'}));
        hold(ax{i},'on');
    end
    %figure('Color','white','Position',[1 -100 1024 768])
    switch lower(plot_type{i})
        case 'energy'
            if ndims(dist.f) == 5
                fprintf('R=%.2f, Phi=%.2f, Z=%.2f\n',eq.fields.r(r0_ind),eq.fields.phi(phi0_ind),eq.fields.z(z0_ind))
                if any(size(dist.f)~=[dist.nenergy dist.npitch dist.nr dist.nz dist.nphi])
                    tmp = squeeze(trapz(dist.phi(phi0_ind:phi1_ind),trapz(dist.z(z0_ind:z1_ind),dist.f,4),5));
                else
                    tmp = squeeze(trapz(dist.phi(phi0_ind:phi1_ind),trapz(dist.z(z0_ind:z1_ind),dist.f(:,:,r0_ind:r1_ind,z0_ind:z1_ind,phi0_ind:phi1_ind),4),5));
                end
            else
                fprintf('R=%.2f,  Z=%.2f\n',eq.fields.r(r0_ind),eq.fields.z(z0_ind))
                if any(size(dist.f)~=[dist.nenergy dist.npitch dist.nr dist.nz dist.nphi])
                    tmp = squeeze(trapz(dist.z(z0_ind:z1_ind),dist.f,4))*2*pi;
                else
                    tmp = squeeze(trapz(dist.z(z0_ind:z1_ind),dist.f(:,:,r0_ind:r1_ind,z0_ind:z1_ind),4))*2*pi;
                end
            end
            rtmp = permute(repmat(dist.r(r0_ind:r1_ind),1,size(tmp,1),size(tmp,2),1),[2,3,1]);
            tmp = squeeze(trapz(dist.r(r0_ind:r1_ind),rtmp.*tmp,3));
            tmp=squeeze(trapz(dist.pitch,tmp,2));
            if fac == 1
                % plot(ax,dist.energy(2:end), tmp(1:end-1),'DisplayName',['Energy - ' name] );
                plot(ax{i},dist.energy, tmp,linestyle,'DisplayName',['Energy - ' name] );
            else
                plot(ax{i},dist.energy, fac*tmp,linestyle,'DisplayName',['Energy - ' name ', scaling factor: ' num2str(fac)]);
            end
            xlabel(ax{i},'Energy [keV]')
            ylabel(ax{i},'Fast Ion Distribution [1/keV]')
        case 'pitch'
            if ndims(dist.f) == 5
                tmp = squeeze(trapz(dphi*nphi/(nphi-1),trapz(dz*nz/(nz-1),trapz(dist.energy,dist.f,1),4),5)); %Integral over energy, z and phi
                tmp = squeeze(trapz(dr*nr/(nr-1),repmat(dist.r(r0_ind:r1_ind)',size(dist.f,2),1).*tmp,2)); %Integral over r with jacobian
            else
                tmp = squeeze(trapz(dz*nz/(nz-1),trapz(dist.energy,dist.f,1),4))*2*pi;
                tmp = squeeze(trapz(dr*nr/(nr-1),repmat(dist.r(r0_ind:r1_ind)',size(dist.f,2),1).*tmp,2));
            end
            if lmean
                tmp=tmp./ sum(vol2d,'all');
            end
            fprintf('Total fast ions in %s: %3.2e\n',file,n_fida);
            plt_data.n_fida=n_fida;
            if fac == 1
                plot(ax{i},dist.pitch, tmp,linestyle,'DisplayName',['Pitch - ' name] );
                fprintf('Total from pitch: %3.2e\n',trapz(dist.pitch, tmp));
            else
                plot(ax{i},dist.pitch, fac*tmp,linestyle,'DisplayName',['Pitch - ' name ', scaling factor: ' num2str(fac)]);
            end
            %plot(dist.pitch, squeeze(trapz(dist.phi,trapz(dist.z,trapz(dist.r,trapz(dist.energy,dist.f,1),3),4),5)),'DisplayName','Pitch');
            xlabel(ax{i},'Pitch [-]')
            ylabel(ax{i},'Fast Ion Distribution [-]')
        case 'bmir'
            disp('bmir is currently not working!')
            continue
            tmp=trapz(dist.energy,dist.f,1);
            modb=sqrt(eq.fields.br.^2+eq.fields.bt.^2+eq.fields.bz.^2);
            i = floor(dist.nr/2);
            j =floor(dist.nz/2);
            k = 1;
            local=dist.f(:,:,i,j,k);
            bmir = modb.*reshape(sign(dist.pitch)./(1-dist.pitch.^2),[1 1 1 numel(dist.pitch)]);
            b=linspace(-5,5,500);
            bmir_ind=discretize(bmir,b);
            bmir_ind(isnan(bmir_ind))=1;
            hist=accumarray(bmir_ind(:),tmp(:));

            bmir=squeeze(bmir(floor(dist.nr/2),floor(dist.nz/2),1,:));
            [bmir, index]=sort(bmir);
            local=local(:,index);
            local = trapz(dist.energy,local,1);
            plot(bmir,local)
            pixplot(dist.energy,bmir,local);
            [maxb,I]=max(modb,[],'all','linear');
            [~,~,index]=ind2sub(size(modb),I);
            minb=min(modb(:,:,index),[],'all');
            yline([-minb minb])
            yline([-modb(floor(dist.nr/2),floor(dist.nz/2),1) modb(floor(dist.nr/2),floor(dist.nz/2),1)])
        case 'epplot'
            if ndims(dist.f) == 5
                tmp = squeeze(trapz(dphi,trapz(dz,dist.f,4),5));
                if ldiff
                    tmp2 = squeeze(trapz(dphi,trapz(dz,dist2.f,4),5));
                end
                rtmp = permute(repmat(dr,1,size(dist.f,1),size(dist.f,2),1),[2,3,1]);
            elseif ndims(dist.f)==4
                tmp = squeeze(trapz(dz,dist.f,4))*2*pi;
                if ldiff
                    tmp2 = squeeze(trapz(dz,dist2.f,4))*2*pi;
                end                
                if r0_ind~=r1_ind
                    rtmp = permute(repmat(rdphi(r0_ind:r1_ind),1,size(dist.f,1),size(dist.f,2),1),[2,3,1]);
                else
                    rtmp = permute(repmat(rdphi(rreq_ind),1,size(dist.f,1),size(dist.f,2),1),[2,3,1]);
                end
            elseif ismatrix(dist.f)
                rtmp = permute(repmat(dr,1,size(dist.f,1),size(dist.f,2),1),[2,3,1]);
                tmp = dist.f;
                if ldiff
                    tmp2 = dist2.f;
                end 
            end
            if r0_ind~= r1_ind
                tmp = squeeze(trapz(dist.r(r0_ind:r1_ind),rtmp.*tmp,3));
                cstring='Fast Ion Distribution [1/keV]';                
                if ldiff
                    tmp2 = squeeze(trapz(dist.r(r0_ind:r1_ind),rtmp.*tmp2,3));
                    if lrel
                        tmp = ((tmp-tmp2)./tmp2).^1;
                        tmp(tmp2<1)=0;
                        cstring='Rel. difference (f_1-f_2)/f_2 [-]'; 
                    else
                        tmp=tmp-tmp2;
                        cstring='Difference (f_1-f_2) [1/keV]'; 
                    end  
                end
            else
                if ndims(tmp)==3
                tmp = squeeze(trapz(dr,rtmp.*tmp,3));
                if ldiff
                    tmp2 = squeeze(trapz(dr,rtmp.*tmp2,3));
                    if lrel
                         tmp = ((tmp-tmp2)./tmp2).^1;
                    else
                        tmp=tmp-tmp2;
                    end  
                end                
                cstring='Fast Ion Distribution [1/keV]';
                else
                    cstring='Local Fast Ion Distribution [1/keV/cm^3]';
                end
            end
            if lmean
                tmp=tmp./sum(vol2d,'all');
            end
            if lcontour
                contour(dist.energy,dist.pitch,tmp',levels,linestyle,'DisplayName',name)
            else
                imagesc(dist.energy,dist.pitch,tmp')
            end
            if lweight
                [~,index(1)]=min(abs(weight.lambda-660));%wvl
                index(2) = channel;
                tmp=squeeze(weight.weight(index(1),:,:,index(2)));
                contour(ax{i},weight.energy,weight.pitch,tmp',levels,linestyle,'DisplayName',name)
            end
            r = dist.r;
            phi = eq.fields.phi;
            z = dist.z;            
            c = colorbar;
            c.Label.String = cstring;
            ylabel('Pitch [-]')
            xlabel('Energy [keV]')
            xlim([dist.energy(1) dist.energy(end)])
            ylim([dist.pitch(1) dist.pitch(end)])
        case 'vflow2d'
            v=sqrt(dist.energy./input.ab/amu*1e3*ec*2);
            tmp=squeeze(trapz(dist.pitch,trapz(v,dist.f.*v,1),2));
            tmp=tmp./dist.denf.*1e-6;
            r = eq.plasma.r;
            z = eq.plasma.z;
            phi=eq.plasma.phi;
            cstring = 'Fast ion vlow vel. [m/s]';
        case 'profiles'
            yyaxis(ax{i},'left')
            plot(ax{i},eq.plasma.r, squeeze(eq.plasma.te(:,zreq_ind,1)), 'DisplayName',['T_e - ' name] );
            hold on
            plot(ax{i},eq.plasma.r, squeeze(eq.plasma.ti(:,zreq_ind,1)), 'DisplayName',['T_i - ' name] );
            plot(ax{i},eq.plasma.r, squeeze(eq.plasma.zeff(:,zreq_ind,1)), 'DisplayName',['Zeff [-] - ' name] );
            plot(ax{i},eq.plasma.r, squeeze(eq.plasma.vt(:,zreq_ind,1)/3e6), 'DisplayName',['Vtor [3e4m/s] - ' name] );
            ylabel(ax{i},'T [keV]')
            legend(ax{i},'Interpreter','none');
            yyaxis(ax{i},'right')
            plot(ax{i},eq.plasma.r, squeeze(eq.plasma.dene(:,zreq_ind,1)), 'DisplayName',['n_e - ' name] );
            xlabel(ax{i},'R [cm]')
            ylabel(ax{i},'n_e [cm^{-3}]')
        case 'profiles_rho'
            yyaxis(ax{i},'left')
            plot(ax{i},eq.plasma.profiles.rho, eq.plasma.profiles.te, 'DisplayName',['T_e - ' name] );
            hold on
            plot(ax{i},eq.plasma.profiles.rho, eq.plasma.profiles.ti, 'DisplayName',['T_i - ' name] );
            plot(ax{i},eq.plasma.profiles.rho, eq.plasma.profiles.zeff, 'DisplayName',['Zeff [-] - ' name] );
            ylabel(ax{i},'T [keV]')
            legend(ax{i},'Interpreter','none');
            yyaxis(ax{i},'right')
            plot(ax{i},eq.plasma.profiles.rho, eq.plasma.profiles.dene, 'DisplayName',['n_e - ' name] );
            xlabel(ax{i},'rho [-]')
            ylabel(ax{i},'n_e [cm^{-3}]')
        case 'ne2d'
            r = eq.plasma.r;
            z = eq.plasma.z;
            phi=eq.plasma.phi;
            tmp = eq.plasma.dene;%.*double(eq.plasma.mask);
            cstring = 'Electron density [m^{-3}]';
        case 'te2d'
            r = eq.plasma.r;
            z = eq.plasma.z;
            phi=eq.plasma.phi;
            tmp = eq.plasma.te;
            cstring = 'Electron temperature [keV]';
        case 'ti2d'
            r = eq.plasma.r;
            z = eq.plasma.z;
            phi=eq.plasma.phi;
            tmp = eq.plasma.ti;
            cstring = 'Ion temperature [keV]';
        case 'vt2d'
            r = eq.plasma.r;
            z = eq.plasma.z;
            phi=eq.plasma.phi;
            tmp = eq.plasma.vt;
            cstring = 'Toroidal Rotation [cm/s]';
        case 'er2d'
            r = eq.plasma.r;
            z = eq.plasma.z;
            phi=eq.plasma.phi;
            tmp = eq.fields.er;
            cstring = 'Electric field in R direcion [V/m]';
        case 'et2d'
            r = eq.plasma.r;
            z = eq.plasma.z;
            phi=eq.plasma.phi;
            tmp = eq.fields.et;
            cstring = 'Electric field in phi direcion [V/m]';
        case 'ez2d'
            r = eq.plasma.r;
            z = eq.plasma.z;
            phi=eq.plasma.phi;
            tmp = eq.fields.ez;
            cstring = 'Electric field in Z direcion [V/m]';
        case 'zeff2d'
            r = eq.plasma.r;
            z = eq.plasma.z;
            phi=eq.plasma.phi;
            tmp = eq.plasma.zeff;
            cstring = 'Effective nuclear charge [-]';
        case 'mask2d'
            r = eq.plasma.r;
            z = eq.plasma.z;
            phi=eq.plasma.phi;
            tmp = double(eq.fields.mask)+double(eq.plasma.mask);
            cstring = 'Boolean mask [-]';
        case 'denn2d'
            r = eq.plasma.r;
            z = eq.plasma.z;
            phi=eq.plasma.phi;
            if isfield(eq.plasma,'denn')
                tmp = eq.plasma.denn;
            else
                disp('No cold neutral density found!')
                continue
            end
            cstring = 'Cold/Edge neutral density [m^{-3}]';
        case 'ba'
            plot(ax{i},eq.plasma.r, squeeze(eq.fields.br(:,zreq_ind,1)),linestyle, 'DisplayName','B_r');
            plot(ax{i},eq.plasma.r, squeeze(eq.fields.bt(:,zreq_ind,1)),linestyle, 'DisplayName','B_t');
            plot(ax{i},eq.plasma.r, squeeze(eq.fields.bz(:,zreq_ind,1)),linestyle, 'DisplayName','B_z');
            xlabel(ax{i},'R [cm]')
            ylabel(ax{i},'Magnetic Field [T]')
            legend(ax{i},'Interpreter','none');
            return;
        case 'fslice'
            disp('fslice is Currently not working!')
            continue
            if index==1
                [~,e_ind]=min(abs(dist.energy-20));
                [~,p_ind]=min(abs(dist.pitch));
                phi_ind=1;
                tmp=squeeze(dist.f(e_ind,p_ind,:,zreq_ind,phi_ind));
            else
                [~,e_ind]=min(abs(dist.energy-index_in(1)));
                [~,p_ind]=min(abs(dist.pitch-index_in(2)));
                [~,zreq_ind]=min(abs(dist.z-index_in(3)));
                [~,phi_ind]=min(abs(dist.phi-index_in(4)));
                tmp=squeeze(dist.f(e_ind,p_ind,:,zreq_ind,phi_ind));
            end
            if fac ==1
                plot(ax{i},dist.r, tmp,linestyle,'DisplayName',['f slice - ' name] );
            else
                plot(ax{i},dist.r, fac*tmp,linestyle,'DisplayName',['f slice - ' name ', scaling factor: ' num2str(fac)]);
            end
            fprintf('Total: %.2e\n', squeeze(trapz(dist.r,tmp)));
            xlabel(ax{i},'R [cm]')
            ylabel(ax{i},'Fast ion distribution slice [1/cm^3/keV/dp]')
            legend(ax{i},'Interpreter','none');
        case 'denf'
            r = dist.r;
            phi = eq.fields.phi;
            z = dist.z;
            tmp = dist.denf;
            if index==1
                index=zreq_ind;
            end
            if fac == 1
                plot(ax{i},dist.r, squeeze(tmp(:,index,1)),linestyle, 'DisplayName',sprintf('%s, z= %.2f',name,dist.z(zreq_ind)));
            else
                plot(ax{i},dist.r, fac*squeeze(tmp(:,index,1)),linestyle, 'DisplayName',['Denf - ' name ', scaling factor: ' num2str(fac)]);
            end
            xlabel(ax{i},'R [m]')
            ylabel(ax{i},'Fast ion density [cm^{-3}]')
            %title(ax{i},sprintf('FI density profile at z= %.2f',dist.z(zreq_ind)))
        case 'fdenf'
            r = dist.r;
            phi = eq.fields.phi;
            z = dist.z;
            if size(dist.f,1)==numel(dist.energy)
                tmp = squeeze(trapz(dist.pitch,trapz(dist.energy,dist.f,1),2));
                if ldiff
                    tmp2 =squeeze(trapz(dist.pitch,trapz(dist.energy,dist2.f,1),2));
                    if lrel
                        tmp = (tmp-tmp2)./tmp2;
                    else
                        tmp=tmp-tmp2;
                    end
                end
            else
                tmp = squeeze(trapz(dist.pitch(p_min:p_max),trapz(dist.energy(e_min:e_max),dist.f(:,:,:,:),1),2));
                if ldiff
                    tmp2 = squeeze(trapz(dist.pitch(p_min:p_max),trapz(dist.energy(e_min:e_max),dist2.f(:,:,:,:),1),2));
                    if lrel
                        tmp = (tmp-tmp2)./tmp;
                    else
                        tmp=tmp-tmp2;
                    end
                end
            end
            if fac == 1
                plot(ax{i},dist.r, tmp(:,zreq_ind),linestyle, 'DisplayName',name );
            else
                plot(ax{i},dist.r, fac*tmp(:,zreq_ind),linestyle, 'DisplayName',[ name ', scaling factor: ' num2str(fac)]);
            end
            disp(['z=',num2str(dist.z(zreq_ind))])
            xlabel(ax{i},'R [m]')
            ylabel(ax{i},'Fast ion density [cm^{-3}]')
            if lrel
            ylabel(ax{i},'Fast ion density (Relative difference)^2 [cm^{-3}]')
            end                
            %title('Fast ion density profile at z=0')
        case 'denf2d'
            r = dist.r;
            z = dist.z;
            phi=dist.phi;
            tmp = dist.denf;%.*dist.r2d;
            if lrel
                tmp=tmp./dist2.denf;
            end
            cstring = 'Fast ion density [cm^{-3}]';
        case {'fdenf2d','fdenftor'}
            r = dist.r;
            phi = dist.phi;
            z = dist.z;
            if size(dist.f,1)~=dist.nenergy
                tmp = squeeze(trapz(dist.pitch(p_min:p_max),trapz(dist.energy(e_min:e_max),dist.f,1),2));
                if lrel
                    tmp = tmp./squeeze(trapz(dist.pitch(p_min:p_max),trapz(dist.energy(e_min:e_max),dist2.f,1),2));
                end
            else
                tmp = squeeze(trapz(dist.pitch(p_min:p_max),trapz(dist.energy(e_min:e_max),dist.f(e_min:e_max,p_min:p_max,:,:,:),1),2));
                if lrel
                    tmp = tmp./squeeze(trapz(dist.pitch(p_min:p_max),trapz(dist.energy(e_min:e_max),dist2.f(e_min:e_max,p_min:p_max,:,:,:),1),2));
                end
            end
            cstring = 'Fast ion density [cm^{-3}]';
        case 'br2d'
            r = eq.fields.r;
            z = eq.fields.z;
            phi = eq.fields.phi;
            tmp = eq.fields.br;
            cstring = 'Magnetic Field B_r [T]';
        case 'bt2d'
            r = eq.fields.r;
            z = eq.fields.z;
            phi = eq.fields.phi;
            tmp = eq.fields.bt;
            cstring = 'Magnetic Field B_t [T]';
        case 'bz2d'
            r = eq.fields.r;
            z = eq.fields.z;
            phi = eq.fields.phi;
            tmp = eq.fields.bz;
            cstring = 'Magnetic Field B_z [T]';
        case 'q2d'
            r = eq.fields.r;
            z = eq.fields.z;
            %Quickly find minor radius
            [~,linearIndMax]=max(eq.plasma.te(:));
            [row,col,~]=ind2sub(size(eq.plasma.te),linearIndMax);
            a2d=sqrt((eq.fields.r2d-eq.fields.r2d(row,col)).^2+(eq.fields.z2d-eq.fields.z2d(row,col)).^2);
            tmp = (eq.fields.bt(:,:,1) .* a2d) ./ sqrt(eq.fields.bz(:,:,1).^2 + eq.fields.br(:,:,1).^2)./eq.fields.r2d;
            cstring = 'Safety factor [-]';
        case 'denftor'
            r = dist.r;
            phi = dist.phi;
            tmp = dist.denf;
            cstring = 'Fast ion density [cm^{-3}]';
            % case 'fdenftor'
            %     r = dist.r;
            %     phi = dist.phi;
            %     z = dist.z;
            %     if numel(index_in)==1
            %         tmp = squeeze(trapz(dist.pitch,trapz(dist.energy,dist.f,1),2));
            %     else
            %         [~,e_min]=min(abs(dist.energy-index_in(1)));
            %         [~,e_max]=min(abs(dist.energy-index_in(2)));
            %         [~,p_min]=min(abs(dist.pitch-index_in(3)));
            %         [~,p_max]=min(abs(dist.pitch-index_in(4)));
            %         tmp = squeeze(trapz(dist.pitch(p_min:p_max),trapz(dist.energy(e_min:e_max),dist.f(e_min:e_max,p_min:p_max,:,:,:),1),2));
            %     end
            %     cstring = 'Fast ion density [m^{-3}]';
        case 'fdenf3d'
            r = dist.r;
            phi = dist.phi;
            z = dist.z;
            if index==1
                tmp = squeeze(trapz(dist.pitch,trapz(dist.energy,dist.f,1),2));
            else
                [~,e_min]=min(abs(dist.energy-index_in(1)));
                [~,e_max]=min(abs(dist.energy-index_in(2)));
                [~,p_min]=min(abs(dist.pitch-index_in(3)));
                [~,p_max]=min(abs(dist.pitch-index_in(4)));
                tmp = squeeze(trapz(dist.pitch(p_min:p_max),trapz(dist.energy(e_min:e_max),dist.f(e_min:e_max,p_min:p_max,:,:,:),1),2));
                phi=[phi; phi(1)];
                tmp=cat(3,tmp,tmp(:,:,1));
            end
            [x_fida,y_fida,z_fida] = ndgrid(r+dr/2,mod(phi,max(phi+dphi/2))+dphi/2,z./1+dz/2);
            denf=permute(tmp(:,:,:),[1 3 2]);
            xg = x_fida .* cos(y_fida);
            yg = x_fida.* sin(y_fida);
            %Patch around is necessary for RPHIZ-XYZ conversion
            p=patch(ax{i},isosurface(xg,yg,z_fida,denf,index));
            colors = parula(20);
            p.FaceColor = colors(10,:);
            p.EdgeColor = 'none';
            view(3)
            axis equal
            camlight
            lighting gouraud
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            cstring = 'Fast ion density [m^{-3}]';
        case 'brtor'
            r = eq.fields.r;
            phi = eq.fields.phi;
            z = eq.fields.z;
            tmp = eq.fields.br;
            cstring = 'Magnetic Field B_r [T]';
        case 'bttor'
            r = eq.fields.r;
            phi = eq.fields.phi;
            z = eq.fields.z;
            tmp = eq.fields.bt;
            cstring = 'Magnetic Field B_t [T]';
        case 'bztor'
            r = eq.fields.r;
            phi = eq.fields.phi;
            z = eq.fields.z;
            tmp = eq.fields.bz;
            cstring = 'Magnetic Field B_z [T]';
        case 'ndens'
            neut_r = sqrt(neut.grid.x_grid(:).^2+neut.grid.y_grid(:).^2);
            [discr,redges]=discretize(neut_r,128);
            r = neut_r;
            tmp = accumarray(discr,neut.nparts(:));
            plot(redges(1:end-1),tmp./diff(redges));
            xlabel('R [m]')
            ylabel('Deposition [particles/s]')
            title('Radial Deposition profile')
        case 'ndens2d'
            %eq.fields.phi=linspace(0,1,32);
            %[r,z,phi] = ndgrid(eq.fields.r,eq.fields.z,eq.fields.phi);
            %uvw = [r(:).*cos(phi(:)), -r(:).*sin(phi(:)),z(:)];
            %uvw = [u(:),v(:),w(:)];
            %[xyz] = uvw_to_xyz(input.alpha, input.beta, input.gamma, uvw, input.origin);
            %ngrid = [neut.grid.x_grid(:),neut.grid.y_grid(:),neut.grid.z_grid(:)];
            neut_r = sqrt(neut.grid.x_grid(:).^2+neut.grid.y_grid(:).^2);
            neut_phi= atan2(neut.grid.y_grid(:),neut.grid.x_grid(:));
            neut_z = neut.grid.z_grid(:);
            
            % [discphi,phiedges]=discretize(neut_phi,32);
            % [discr,redges]=discretize(neut_r,64);
            % [discz,zedges]=discretize(neut_z,128);
            % r = redges(1:end-1)+mean(diff(redges))/2;
            % phi = phiedges(1:end-1)+mean(diff(phiedges))/2;
            % z = zedges(1:end-1)+mean(diff(zedges))/2;
            r=linspace(min(neut_r,[],'all'),max(neut_r,[],'all'),64);
            phi=linspace(min(neut_phi,[],'all'),max(neut_phi,[],'all'),32);
            z=linspace(min(neut_z,[],'all'),max(neut_z,[],'all'),64);
            %plot2DHistogram(neut.dens(:),neut_r,neut_z,[],[],'nres',32);
            %tmp = accumarray([discr,discz],neut_r.*neut.dens(:));
            [rg,zg,phig]=ndgrid(r,z,phi);
            dens=scatteredInterpolant(neut_r,neut_z,neut_phi,neut.dens(:),'linear','none');
            tmp=dens(rg,zg,phig);
            tmp=sum(rg.*tmp,3,'omitmissing')*(phi(2)-phi(1));%Sum over phi
            %             ndens_F = scatteredInterpolant(ngrid,neut.dens(:));
            %             ndens_F.Method = 'linear';
            %             ndens_F.ExtrapolationMethod = 'none';
            %             ndens = ndens_F(uvw);
            %             tmp=reshape(ndens,size(r));
            cstring='Neutral Density [neutrals/cm^3]';
            index=1;
        case 'ndenstor'
            neut_r = sqrt(neut.grid.x_grid(:).^2+neut.grid.y_grid(:).^2);
            neut_phi= atan2(neut.grid.y_grid(:),neut.grid.x_grid(:));
            neut_z = neut.grid.z_grid(:);
            [discphi,phiedges]=discretize(neut_phi,eq.fields.nphi);
            [discr,redges]=discretize(neut_r,eq.fields.nr);
            [discz,zedges]=discretize(neut_z,eq.fields.nz);
            r = redges(1:end-1)+mean(diff(redges))/2;
            phi = phiedges(1:end-1)+mean(diff(phiedges))/2;
            z = zedges(1:end-1)+mean(diff(zedges))/2;
            tmp = accumarray([discr,discz,discphi],neut_r.*neut.dens(:));
            %tmp=sum(tmp,3)*(phiedges(2)-phiedges(1));%Sum over phi
            %             ndens_F = scatteredInterpolant(ngrid,neut.dens(:));
            %             ndens_F.Method = 'linear';
            %             ndens_F.ExtrapolationMethod = 'none';
            %             ndens = ndens_F(uvw);
            %             tmp=reshape(ndens,size(r));
            cstring='Neutral Density [neutrals/cm^3]';
        case 'ndensvert'
            pixplot(neut.grid.x, neut.grid.z, squeeze(sum(neut.tdens(:,:,index,:) + neut.hdens(:,:,index,:) + neut.fdens(:,:,index,:), 1)));
            xlabel('Beam Grid X [cm]')
            ylabel('Beam Grid Z [cm]')
            cstring = 'Beam neutral density [1/cm^3]';
            c = colorbar;
            c.Label.String = cstring;
        case 'ndenshorz'
            pixplot(neut.grid.x, neut.grid.y, squeeze(sum(neut.tdens(:,:,:,index) + neut.hdens(:,:,:,index)+ neut.fdens(:,:,:,index), 1)));
            xlabel('Beam Grid X [cm]')
            ylabel('Beam Grid Y [cm]')
            cstring = 'Beam neutral density [1/cm^3]';
            c = colorbar;
            c.Label.String = cstring;
        case 'ndenscross'
            pixplot(neut.grid.y, neut.grid.z, squeeze(sum(neut.tdens(:,index,:,:) + neut.hdens(:,index,:,:)+ neut.fdens(:,index,:,:), 1)));
            xlabel('Beam Grid Y [cm]')
            ylabel('Beam Grid Z [cm]')
            cstring = 'Beam neutral density [1/cm^3]';
            c = colorbar;
            c.Label.String = cstring;
        case 'fdensvert'
            pixplot(neut.grid.x, neut.grid.z, squeeze(sum( neut.fdens(:,:,index,:), 1)));
            xlabel('Beam Grid X [cm]')
            ylabel('Beam Grid Z [cm]')
            cstring = 'Beam neutral density [1/cm^3]';
            c = colorbar;
            c.Label.String = cstring;
        case 'fdenshorz'
            pixplot(neut.grid.x, neut.grid.y, squeeze(sum(neut.fdens(:,:,:,index), 1)));
            xlabel('Beam Grid X [cm]')
            ylabel('Beam Grid Y [cm]')
            cstring = 'Beam neutral density [1/cm^3]';
            c = colorbar;
            c.Label.String = cstring;
        case 'fdenscross'
            pixplot(neut.grid.y, neut.grid.z, squeeze(sum(neut.fdens(:,index,:,:), 1)));
            xlabel('Beam Grid Y [cm]')
            ylabel('Beam Grid Z [cm]')
            cstring = 'Beam neutral density [1/cm^3]';
            c = colorbar;
            c.Label.String = cstring;            
        case 'halovert'
            pixplot(neut.grid.x, neut.grid.z, squeeze(sum(neut.halodens(:,:,index,:) + neut.dcxdens(:,:,index,:), 1)));
            xlabel('Beam Grid X [cm]')
            ylabel('Beam Grid Z [cm]')
            cstring = 'Halo+DCX neutral density [1/cm^3]';
            c = colorbar;
            c.Label.String = cstring;
            set(ax{i},'ColorScale','log')
        case 'halohorz'
            pixplot(neut.grid.x, neut.grid.y, squeeze(sum(neut.halodens(:,:,:,index) + neut.dcxdens(:,:,:,index), 1)));
            xlabel('Beam Grid X [cm]')
            ylabel('Beam Grid Y [cm]')
            cstring = 'Halo+DCX neutral density [1/cm^3]';
            c = colorbar;
            c.Label.String = cstring;
            set(ax{i},'ColorScale','log')
        case 'halocross'
            pixplot(neut.grid.y, neut.grid.z, squeeze(sum(neut.halodens(:,index,:,:) + neut.dcxdens(:,index,:,:), 1)));
            xlabel('Beam Grid Y [cm]')
            ylabel('Beam Grid Z [cm]')
            cstring = 'Halo+DCX neutral density [1/cm^3]';
            c = colorbar;
            c.Label.String = cstring;
            set(ax{i},'ColorScale','log')
        case {'weights','weights_dist','weight_dist'}
            [~,index_in(1)]=min(abs(weight.lambda-index_in(1)));%wvl
            if  numel(index_in)==1
                index_in(2) = 1;%channel
            end
            tmp=squeeze(weight.weight(index_in(1),:,:,index_in(2)));
            if strcmp(plot_type{i}(end-3:end),'dist')
                tmp2=squeeze(weight.mean_f(:,:,index_in(2)));
            end
            if lcontour
                contour(ax{i},weight.energy,weight.pitch,tmp',levels,linestyle,'DisplayName',name)
            else
                imagesc(ax{i},weight.energy,weight.pitch,tmp');

            end
            if ~isempty(tmp2)
                currentLimits=clim;
                contour(ax{i},weight.energy,weight.pitch,tmp2',5,linestyle,'LineColor','w','DisplayName','FI Dist.')
                clim(currentLimits);
            end
            c = colorbar(ax{i});
            c.Label.String = 'Sensitivity [ph*cm/(s*fast ion)]';
            xlabel('Energy [keV')
            ylabel('Pitch [-]')
            ylim([-1 1])
            title(sprintf('%s, R=%.2fcm',char(deblank(geom.spec.id(index_in(2)))), geom.spec.radius(index_in(2))));
        case 'spectrum'
            specr = spec.full + spec.half + spec.third + spec.halo + spec.dcx + spec.fida;% + spec.brems;
            if isfield(spec,'pfida')
                specr=specr+spec.pfida;
            else
                disp('No passive FIDA signal found!')
            end
            specr = specr.*fac;
            if lmean
                k = 12;
                specr = movmean(specr,k);
                disp(['Applying moving mean with length ', num2str(k), ' to FIDASIM data: ', file]);
            end
            if ~isempty(sim_data)
                %                 for j = 1:size(sim_data.lambda,2)
                %                     spectmp(:,j) = interp1(spec.lambda, specr(:,j), sim_data.lambda(:,j),'pchip');
                %                 end
                %                 disp('Interpolated wavelength to match data');
                if fac~=1.0
                    name = [name,', scale=',num2str(fac)];
                end
                %cwav_mid=sim_data.cwav_mid(channel);
                cwav_mid=mean(spec.lambda);%-(spec.lambda(2)-spec.lambda(1));
                %cwav_mid = interp1(1:size(spec.lambda,1),spec.lambda,size(spec.lambda,1)/2.)-(spec.lambda(2)-spec.lambda(1));

                disp(['Cwav_mid_fidasim=', num2str(cwav_mid)]);
                %cwav_mid=sim_data.cwav_mid(channel);
                %cwav_mid = interp1(1:size(spec.lambda,1),spec.lambda,size(spec.lambda,1)/2.)-(spec.lambda(2)-spec.lambda(1))/2.;
                % We need to flip the kernel to do the same thing as
                % fplot...
                instfu = flipud(box_gauss_funct(spec.lambda,0.,1.,cwav_mid,sim_data.instfu_gamma,sim_data.instfu_box_nm));

                plot(spec.lambda,conv(specr(:,channel),instfu(:,channel),'same'), 'DisplayName', ['Spectrum - ' name] );
                %plot(spec.lambda,specr(:,channel), 'DisplayName', ['Spectrum no instfu - ' name] );
                %plot(spec.lambda,specr(:,channel), 'DisplayName', ['Spectrum - ' name] );
                %plot(spec.lambda, conv(spec.full(:,channel),instfu(:,channel),'same'), 'DisplayName',['Full - ' name] );
                %plot(spec.lambda, spec.full(:,channel), 'DisplayName',['Full - ' name] );
                %plot(spec.lambda, conv(spec.half(:,channel),instfu(:,channel),'same'),  'DisplayName',['Half - ' name] );
                %plot(spec.lambda, conv(spec.third(:,channel),instfu(:,channel),'same'),  'DisplayName',['Third - ' name] );
                if ( isfield(spec,'pfida') && lpassive)
                    plot(spec.lambda, conv(spec.pfida(:,channel),instfu(:,channel),'same'),  'DisplayName',['Passive FIDA - ' name] );
                end
                if ( isfield(spec,'brems') && lbrems)
                    plot(spec.lambda, conv(spec.brems(:,channel),instfu(:,channel),'same'),  'DisplayName',['Bremsstrahlung - ' name] );
                end
                %plot(spec.lambda, conv(spec.halo(:,channel)+spec.dcx(:,channel),instfu(:,channel),'same'),  'DisplayName',['Halo+DCX - ' name] ); %+spec.brems(:,channel)
                %plot(spec.lambda, conv(spec.dcx(:,channel),instfu(:,channel),'same'),  'DisplayName',['DCX only - ' name] ); %+spec.brems(:,channel)
                %fprintf('Halo Centered at %.3f nm\n', sum(spec.lambda.*conv(spec.halo(:,channel)+spec.dcx(:,channel),instfu(:,channel),'same'))./sum(conv(spec.halo(:,channel)+spec.dcx(:,channel),instfu(:,channel),'same')));
                %plot(spec.lambda, conv(spec.fida(:,channel),instfu(:,channel),'same'),  'DisplayName',['FIDA - ' name] );
            else
                plot(spec.lambda,specr(:,channel),linestyle, 'DisplayName', ['Spectrum - ' name] );
                if (isfield(spec,'pfida' ) && lpassive)
                    plot(spec.lambda, spec.pfida(:,channel),  'DisplayName',['Passive FIDA - ' name] );
                end
                if ( isfield(spec,'brems') && lbrems)
                    plot(spec.lambda,spec.brems(:,channel),  'DisplayName',['Bremsstrahlung - ' name] );
                end
                plot(spec.lambda,spec.fida(:,channel),linestyle, 'DisplayName', ['FIDA - ' name] );
                %plot(spec.lambda, spec.halo(:,channel)+spec.dcx(:,channel),  'DisplayName',['Halo+DCX - ' name] ); %+spec.brems(:,channel)
                plot(spec.lambda, spec.halo(:,channel),  'DisplayName',['Halo - ' name] ); %+spec.brems(:,channel)
                plot(spec.lambda, spec.dcx(:,channel),  'DisplayName',['DCX - ' name] ); %+spec.brems(:,channel)
                %plot(spec.lambda, spec.full(:,channel),  'DisplayName',['Full - ' name] ); %+spec.brems(:,channel)
                %plot(spec.lambda, spec.half(:,channel),  'DisplayName',['Half - ' name] ); %+spec.brems(:,channel)
                %plot(spec.lambda, spec.third(:,channel),  'DisplayName',['Third - ' name] ); %+spec.brems(:,channel)
                disp('Supply in_data from e.g. get_bes_fida_aug_data for more plots!')
            end
            hold on
            xlabel('Wavelength [nm]')
            ylabel('Intensity [Ph/(s nm m^2 sr)]')
            set(gca,'YScale','log')
            xlim([650 663]);
            %ylim([1e15, 3e19]);
            if lgeom
                %title(['Channel: ' geom.spec.id(channel)])
                disp(['Channel: ', char(geom.spec.id(channel))])
                disp(['R= ', num2str(geom.spec.radius(channel))])
            end
            %legend(ax{i},'Interpreter','none','Location','northeast');
        case 'los3d'
            vec = [0, 0, -1];
            lens = rotate_points(geom.spec.lens,vec,deg2rad(rotation));
            axi = rotate_points(geom.spec.axis,vec,deg2rad(rotation));
            %los = [lens, lens + axi.*max(sqrt(sum(lens(channel,1:2).^2,2)))*length];
            %los = reshape(los,geom.spec.nchan,3,2);
            [xpts,ypts,zpts] = getPointsAlongAxis(lens', axi', 50, length);
            plot3(ax{i},xpts(:,channel)*fac,ypts(:,channel)*fac,zpts(:,channel)*fac,linestyle);
            src = rotate_points(geom.nbi.src,vec,deg2rad(rotation));
            axi_nbi = rotate_points(geom.nbi.axis,vec,deg2rad(rotation));
            [xpts,ypts,zpts] = getPointsAlongAxis(src', axi_nbi', 50, length);
            plot3(ax{i},xpts*fac,ypts*fac,zpts*fac,'-r');          
            plot3(ax{i},xpts(1)*fac,ypts(1)*fac,zpts(1)*fac,'+k');  
            %los_nbi = [src, src + axi_nbi.*length.*sqrt(sum(src(1:2).^2))];
            % los_nbi = reshape(los_nbi,3,2);
            % plot3(ax{i},squeeze(los_nbi(1,:))'*fac,squeeze(los_nbi(2,:))'*fac,squeeze(los_nbi(3,:))'*fac,'-r');
            % plot3(ax{i},squeeze(los_nbi(1,1))'*fac,squeeze(los_nbi(2,1))'*fac,squeeze(los_nbi(3,1))'*fac,'+k');
            if isfield(input,'xmin')
                coords = [input.xmin input.ymin input.zmin;...
                    input.xmax input.ymin input.zmin;...
                    input.xmax input.ymax input.zmin;...
                    input.xmin input.ymax input.zmin;...
                    input.xmin input.ymin input.zmax;...
                    input.xmax input.ymin input.zmax;...
                    input.xmax input.ymax input.zmax;...
                    input.xmin input.ymax input.zmax];
                %UVW Corrdinates are Machine coordinates
                coords= xyz_to_uvw(input.alpha, input.beta, input.gamma, coords, input.origin);
                faces = [1 2 3 4 1;
                    1 2 6 5 1;
                    2 3 7 6 2;
                    3 4 8 7 3;
                    4 1 5 8 4;
                    5 6 7 8 5];
                hold on;
                ha = patch('vertices',coords,'faces',faces);
                set(ha,'FaceColor','red', 'facealpha', 0.1);
                plot3(coords(:,1),coords(:,2),coords(:,3),'.','DisplayName','Edge Points')
                plot3(input.origin(1),input.origin(2),input.origin(3),'+','DisplayName','Beam Grid Origin')
            end
            if isfield(eq,'fields')
                if eq.fields.nphi==1
                    eq.fields.phi=0:0.2:2*pi;
                end
                [r,phi,z]=ndgrid([eq.fields.r(1),eq.fields.r(end)],eq.fields.phi,[eq.fields.z(1),eq.fields.z(end)]);
                x=r.*cos(phi);
                y=r.*sin(phi);
                k=boundary(x(:),y(:),z(:));
                trisurf(k,x,y,z,'Facecolor','red','EdgeColor','none','FaceAlpha',0.1)
                if lsep
                    dphi=eq.fields.phi(2)-eq.fields.phi(1);
                    dr=(eq.fields.r(2)-eq.fields.r(1))*fac;
                    dz=(eq.fields.z(2)-eq.fields.z(1))*fac;
                    [r,phi,z_fida] = ndgrid(eq.fields.r*fac+dr/2,eq.fields.phi+dphi/2,eq.fields.z*fac+dz/2);
                    x_fida=r.*cos(phi);
                    y_fida=r.*sin(phi);
                    tmp=permute(eq.plasma.dene,[1 3 2]);
                    if eq.fields.nphi==1
                        tmp=repmat(tmp,1,numel(eq.fields.phi),1);
                    end
                    N=scatteredInterpolant(x_fida(:),y_fida(:),z_fida(:),tmp(:),'linear','none');
                    xg=linspace(min(x_fida,[],'all'),max(x_fida,[],'all'),51);
                    yg=linspace(min(y_fida,[],'all'),max(y_fida,[],'all'),52);
                    zg=linspace(min(z_fida,[],'all'),max(z_fida,[],'all'),53);
                    [x_dist,y_dist,z_dist] = meshgrid(xg,yg,zg);
                    dene=N(x_dist,y_dist,z_dist);
                    isosurface(x_dist,y_dist,z_dist,dene,4e13);
                end
            end
            xlabel('X [cm]')
            ylabel('Y [cm]')
            zlabel('Z [cm]')
            camlight left
            %sname = [filename, '_', plot_type{i}];
            %writematrix(los_nbi,sname,linestyle);
            % set(h, {'DisplayName'}, cellstr(deblank(geom.spec.id(channel))))
            %legend(h,'Location','bestoutside');
            %axis equal;
            % rotate3d on;
        case 'lostor'
            vec = [0, 0, -1];
            lens = rotate_points(geom.spec.lens',vec,deg2rad(rotation))'; %AUG: 67.5
            axi = rotate_points((geom.spec.lens + geom.spec.axis.*max(geom.spec.radius)*length)',vec,deg2rad(rotation))';
            los = [lens, axi];
            los = reshape(los,3,geom.spec.nchan,2);
            for k = 1:size(channel,2)
                %,'Color',ax{i}.ColorOrder(k,:)
                if iscell(chan_description)
                    displ=chan_description{k};
                else
                    displ=chan_description;
                end
                h=plot(ax{i},squeeze(los(1,channel(:,k),:))'*fac,squeeze(los(2,channel(:,k),:))'*fac ,linestyle, 'DisplayName', displ);
                %set(h, 'DisplayName', chan_description{k});
            end
            %plot(ax{i},squeeze(geom.spec.closest_points(1,:))'*fac,squeeze(geom.spec.closest_points_cyl(2,:))'*fac,'dk');
            %h=plot(ax,squeeze(los(1,channel,:))'*fac,squeeze(los(2,channel,:))'*fac, 'k');
            %set(h, {'DisplayName'}, cellstr(deblank(geom.spec.id(channel))));
            %legend(h,'Location','bestoutside');
            xlabel('X [cm]')
            ylabel('Y [cm]')
            return
        case 'los2d'
            vec = [0, 0, -1];
            %lens = rotate_points(geom.spec.lens',vec,deg2rad(67.5))';
            %axi = rotate_points(geom.spec.axis',vec,deg2rad(67.5))';
            lens = rotate_points(geom.spec.lens',vec,deg2rad(rotation))';
            axi = rotate_points(geom.spec.axis',vec,deg2rad(rotation))';
            los = [lens, lens + axi.*max(geom.spec.radius)*length];
            los = reshape(los,3,geom.spec.nchan,2);
            lost = permute(los,[3,2,1]);
            los2 = reshape(lost,2,geom.spec.nchan*3);
            losre=interp1([0,1],los2,linspace(0,1,200));
            los= reshape(losre,[],size(lost,2),size(lost,3));
            for k = 1:size(channel,2)
                r = sqrt(los(:,channel(:,k),1).^2 + los(:,channel(:,k),2).^2);
                %,'Color',ax{i}.ColorOrder(k,:)
                if iscell(chan_description)
                    displ=chan_description{k};
                else
                    displ=chan_description;
                end
                displ=geom.spec.id(channel);
                displ=deblank(displ(k));
                h=plot(ax{i},r*fac,squeeze(los(:,channel(:,k),3))*fac,linestyle, 'DisplayName', displ);
                %set(h, 'DisplayName', chan_description{k});
            end
            return
            %legend(h,'Location','bestoutside');
        case {'birth_r','birth_r_gc'}
            if strcmp(plot_type{i}(end-1:end),'gc')
                edges = min(birth.ri_gc(1,:)):1:max(birth.ri_gc(1,:));
                dists = discretize(birth.ri_gc(1,:),edges);
            else
                edges = min(birth.ri(1,:)):1:max(birth.ri(1,:));
                dists = discretize(birth.ri(1,:),edges);
            end
            dists(isnan(dists)) = 1;
            weights = birth.weight;
            %sum(weights)
            histo = accumarray(dists',weights,[size(edges,2)-1, 1]);
            x=(edges(2:end-1)+mean(diff(edges))/2)/100;
            plot(ax{i},x,histo(2:end),'DisplayName','FIDASIM','LineWidth', 2.0)
            xlabel('R [m]')
            ylabel('Deposition [particles/s]')
        case {'birth_z','birth_z_gc'}
            if strcmp(plot_type{i}(end-1:end),'gc')
                edges = min(birth.ri_gc(2,:)):1:max(birth.ri_gc(2,:));
                dists = discretize(birth.ri_gc(2,:),edges);
            else
                edges = min(birth.ri(2,:)):1:max(birth.ri(2,:));
                dists = discretize(birth.ri(2,:),edges);
            end
            dists(isnan(dists)) = 1;
            weights = birth.weight;
            %sum(weights);
            histo = accumarray(dists',weights,[size(edges,2)-1, 1]);
            x=(edges(2:end-1)+mean(diff(edges))/2)/100;
            plot(ax{i},x,histo(2:end),'DisplayName','FIDASIM','LineWidth', 2.0)
            xlabel('Z [m]')
            ylabel('Deposition [particles/s]')
        case {'birth_phi','birth_phi_gc'}
            if strcmp(plot_type{i}(end-1:end),'gc')
                birth_phi=mod(birth.ri_gc(3,:),2*pi);
                edges = min(birth_phi):0.01:max(birth_phi);
                dists = discretize(birth_phi,edges);
            else
                birth_phi=mod(birth.ri(3,:),2*pi);
                edges = min(birth_phi):0.01:max(birth_phi);
                dists = discretize(birth_phi,edges);
            end
            dists(isnan(dists)) = 1;
            weights = birth.weight;
            %sum(weights)
            histo = accumarray(dists',weights,[size(edges,2)-1, 1]);
            x=(edges(2:end-1)+mean(diff(edges))/2);
            plot(ax{i},x,histo(2:end),'DisplayName','FIDASIM','LineWidth', 2.0)
            xlabel('Phi [rad]')
            ylabel('Deposition [particles/s]')
        case 'birth_pitch'
            edges=linspace(-1.05,1.05,70)';
            dists = discretize(birth.pitch,edges);
            dists(isnan(dists)) = 1;
            weights = birth.weight;
            sum(weights);
            histo = accumarray(dists,weights',[numel(edges)-1, 1]);
            x=(edges(2:end-1)+mean(diff(edges))/2);
            plot(ax{i},x,histo(2:end),'DisplayName','FIDASIM','LineWidth', 2.0)
            xlabel('Phi [rad]')
            ylabel('Deposition [particles/s]')
    end
    %disp(plot_type{i});
    %if numel(plot_type{i}) > 2
    if strcmp(plot_type{i}(end-1:end),'2d')
        if ltorint
            if nphi>1
                tmp = trapz(dphi*nphi/(nphi-1),dist.r2d.*tmp,3);
            else
                tmp=dist.r2d.*tmp*2*pi;
            end
            index=1;
            cstring(end-2)='2';%Denote area density
        end
        if size(tmp,3)==1&&ldist
            index=1;
        end
        if lcontour
            contour(ax{i},r*fac,z*fac,squeeze(tmp(:,:,index))',levels,linestyle,'DisplayName',name)
        else
            imagesc(ax{i},r*fac,z*fac,tmp(:,:,index)');
            c = colorbar(ax{i});
            c.Label.String = cstring;
        end
        if ~isempty(dist) &&~ltorint
            if dist.nphi > 1
                title(ax{i},sprintf('phi=%.2f',dist.phi(phireq_ind)))
            end
        end
        if lsep
            contour(ax{i},eq.plasma.r*fac,eq.plasma.z*fac,squeeze(eq.plasma.dene(:,:,phireq_ind))',[1 1],'w-','DisplayName','')
        end
        if lintersection
            if channel==0
                intersections = calculateIntersections(geom.spec.lens, geom.spec.axis, phi(phireq_ind));
            else
                intersections = calculateIntersections(geom.spec.lens(:,channel), geom.spec.axis(:,channel), phi(phireq_ind));
            end
            plot(ax{i},intersections(1,:),intersections(2,:),'k.');
        end
        if strcmp(get(ax{i},'XLimMode'),'auto')
            xlim(ax{i},[r(1)*fac r(end)*fac])
            ylim(ax{i},[z(1)*fac z(end)*fac])
            axis equal
        end
        if fac==1
            xlabel(ax{i},'R [cm] ')
            ylabel(ax{i},'Z [cm] ')
        else
            xlabel(ax{i},['R [cm] * ', num2str(fac)])
            ylabel(ax{i},['Z [cm] * ', num2str(fac)])
        end

    elseif strcmp(plot_type{i}(end-2:end),'tor') && ldist
        if index==1
            index=zreq_ind;
        end
        if ndims(dist.f) < 5
            disp('4D Distribution has no toroidal information')
            return;
        end
        if lsep
            contour(ax{i},eq.plasma.r*fac,eq.plasma.phi*fac,squeeze(eq.plasma.dene(:,zreq_ind,:))',[1e11 1e11],'w-','DisplayName','');
        end
        % Shift theta values
        phi_shifted = mod(phi + pi, 2*pi) - pi;

        % Sort the data according to the shifted theta
        [phi, idx] = sort(phi_shifted);
        if ~ismatrix(tmp)
            tmp_shifted = squeeze(tmp(:,zreq_ind,:));
        else
            tmp_shifted = squeeze(tmp);
        end
        tmp_shifted=tmp_shifted(:, idx);
        % imagesc(r,phi,tmp_shifted');
        if lcontour
            contour(ax{i},r*fac,phi,tmp_shifted',levels,linestyle,'DisplayName',name);
        else
            % [tmpx,tmpy]=meshgrid(r*fac,phi);
            % [X,Y]=pol2cart(tmpy,tmpx);
            % % Create the interpolant
            % F = scatteredInterpolant(X(:), Y(:), tmp_shifted(:), 'linear', 'none');
            %
            % % Define the grid for interpolation
            % xq = linspace(min(X(:)), max(X(:)), 200);
            % yq = linspace(min(Y(:)), max(Y(:)), 200);
            % [Xq, Yq] = meshgrid(xq, yq);
            %
            % % Interpolate the data onto the grid
            % %tmp_shifted = F(Xq, Yq);
            % % Define the vertices
            % vertices = [X(:), Y(:)];
            %
            % % Define the faces (connectivity)
            % faces = [];
            % for k = 1:size(tmpx, 1) - 1
            %     for j = 1:size(tmpx, 2) - 1
            %         v1 = (k-1)*size(tmpx, 2) + j;
            %         v2 = v1 + 1;
            %         v3 = v1 + size(tmpx, 2) + 1;
            %         v4 = v1 + size(tmpx, 2);
            %         faces = [faces; v1, v2, v3, v4];
            %     end
            % end
            % patch('Faces', faces, 'Vertices', vertices, 'FaceVertexCData', tmp_shifted(:), 'FaceColor', 'interp', 'EdgeColor', 'none');
            imagesc(ax{i},r*fac,phi,tmp_shifted');
            c = colorbar(ax{i});
            c.Label.String = cstring;
        end
        if lintersection
            if channel==0
                for j=1:size(geom.spec.lens,2)
                    [intersections(1,:,j),intersections(2,:,j),intersections(3,:,j)] = getPointsAlongAxis(geom.spec.lens(:,j), geom.spec.axis(:,j),50,500,true);
                end
            else
                chandex=find(channel);
                for j=1:numel(chandex)
                    [intersections(1,:,j),intersections(2,:,j),intersections(3,:,j)] = getPointsAlongAxis(geom.spec.lens(:,chandex(j)), geom.spec.axis(:,chandex(j)),50,500,true);
                end
            end
            [intersections(1,:,j+1),intersections(2,:,j+1),intersections(3,:,j+1)]=getPointsAlongAxis(geom.nbi.src, geom.nbi.axis,50,1000,true);
            % r_coords=squeeze(sqrt(sum(intersections(1:2,:,:).^2,1)));
            % phi_coords=squeeze(atan2(intersections(2,:,:),intersections(1,:,:)));
            %plot(ax{i},r_coords,phi_coords);
            plot(ax{i},squeeze(intersections(1,:,:)),squeeze(intersections(2,:,:)));
        end

        xlabel(ax{i},'R [cm]')
        ylabel(ax{i},'Phi [rad]')
        title(ax{i},sprintf('Z=%.2fcm',dist.z(zreq_ind)))
        xlim(ax{i},[r(1) r(end)])
        ylim(ax{i},[phi(1) phi(end)])
    elseif strcmp(plot_type{i}(end-2:end),'tor')
        if lcontour
            contour(ax{i},r*fac,phi,squeeze(tmp(:,zreq_ind,:))',levels,linestyle,'DisplayName',name);
        else
            imagesc(ax{i},r*fac,phi,squeeze(tmp(:,zreq_ind,:))');
            c = colorbar(ax{i});
            c.Label.String = cstring;
        end
        %pixplot(r,phi,squeeze(tmp(:,index,:)))
        xticks(unique(round(r,2,'significant')))
        yticks(unique(round(phi,2,'significant')))
        xlabel('R [cm]')
        ylabel('Phi [rad]')
        title(ax{i},sprintf('Z=%.2fcm',z(zreq_ind)))

        xlim([r(1) r(end)])
    end

    if liota
        [src(1),src(2),src(3)]=pol2cart(phi(phireq_ind),r(end),z(zreq_ind));
        [ptax(1),ptax(2),ptax(3)]=pol2cart(phi(phireq_ind),r(end)-1,z(zreq_ind));
        ptax=ptax-src;
        [pts(1,:),pts(2,:),pts(3,:)]=getPointsAlongAxis(src, ptax,800,r(end)-r(1),true);
        x_plt={};
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
            qiota='q';
        elseif ~isempty(vmec)
            if ~isfield(vmec,'Fchi')
                vmec=vmec_rzphi_s_interp(vmec);
            end
            phin = vmec.phi./vmec.phi(end);
            s_line = vmec.Fchi(pts(1,:)/100,mod(pts(2,:),vmec.zeta(end)),pts(3,:)/100);
            q_line=interp1(phin,vmec.iotaf,s_line);
            for m=1:numel(iota_vec)
                [~,tmpp]=findpeaks(-abs(q_line-iota_vec(m)),'MinPeakHeight',-.01);
                if ~isempty(tmpp)
                    x_plt{m}=tmpp;
                else
                    disp(['Found no surfaces for q=', num2str(iota_vec(m))])
                end
            end
            qiota='\iota';
        end
        % if lrho
        %     if ~isempty(vmec)
        %         pts(1,:) = sqrt(vmec.Fchi(pts(1,:)/100,mod(pts(2,:),vmec.zeta(end)),pts(3,:)/100));
        %     else
        %         pts(1,:) = interp3(fida_data.eq.fields.r,fida_data.eq.fields.phi,fida_data.eq.fields.z,...
        %             permute(sqrt(fida_data.eq.fields.s),[ 3 1 2]),pts(1,:),mod(pts(2,:),fida_data.eq.fields.phi(end)),pts(3,:),'linear',NaN);
        %     end
        % end
        for m=1:numel(iota_vec)
            xl=xline(ax{i},pts(1,x_plt{m}),'-',{sprintf('%s=%.2f',qiota,iota_vec(m))},'HandleVisibility','off');
            for j=1:numel(xl)
                xl(j).LabelHorizontalAlignment='center';
                xl(j).LabelVerticalAlignment='bottom';
            end
        end
    end

    if lsave
        if llegend
            legend(ax{i},'Interpreter','none');
        end
        sname = [file, '_', name,  '_', plot_type{i},'_',sprintf('%d',index_in)];
        savefig(ax{i}.Parent,[sname,'.fig'])
        set(ax{i}.Parent, 'Renderer', 'painters');
        set(ax{i}, 'Color', 'none');
        if leps
            exportgraphics(ax{i}.Parent,[sname,'.eps'],'Resolution',300,'BackgroundColor','none');
        end
        exportgraphics(ax{i}.Parent,[sname,'.png'],'Resolution',600);
    end
    plt_data.(plot_type{i})=tmp;
    plt_data.zreq_ind=zreq_ind;
    if zreq_ind~=0; plt_data.zloc=z(zreq_ind); end
    plt_data.phireq_ind=phireq_ind;
    if phireq_ind ~=0; plt_data.philoc=phi(phireq_ind); end
    plt_data.rreq_ind=rreq_ind;
    if rreq_ind~=0; plt_data.rloc=r(rreq_ind); end
    plt_data.index=index;
    
end


end


function [dist, phireq_ind, eq, z0_ind, zreq_ind, phi1_ind, zval, r0_ind, r1_ind, z1_ind, phi0_ind] = load_dist_and_ind(index_in, dist, dist_name, phireq_ind, eq, rreq_ind, z0_ind, zreq_ind, e_min, p_min, e_max, p_max, ltor, phi1_ind, lz, zval, ltrim, r0_ind, r1_ind, z1_ind, phi0_ind)
if isempty(index_in) %If not specified, use whole distribution!
    dist.f= h5read(dist_name,'/f');
    dist.denf= h5read(dist_name,'/denf');
elseif numel(index_in)==1
    phireq_ind=index_in;%
    fprintf('Phi=%.2f\n',eq.fields.phi(phireq_ind))
    dist.f= h5read(dist_name,'/f',[1,1,1,1,phireq_ind],[Inf Inf Inf Inf 1]);
    dist.denf= h5read(dist_name,'/denf',[1,1,phireq_ind],[Inf Inf 1]);
elseif numel(index_in)==3
    fprintf('R=%.2f, Phi=%.2f, Z=%.2f\n',eq.fields.r(rreq_ind),eq.fields.phi(phireq_ind),eq.fields.z(zreq_ind))
    dist.f= h5read(dist_name,'/f',[1,1,rreq_ind,zreq_ind,phireq_ind],[Inf Inf 1 1 1]);
    dist.denf= h5read(dist_name,'/denf',[rreq_ind,zreq_ind,phireq_ind], [1 1 1]);
    r0_ind=rreq_ind;
    r1_ind=rreq_ind;
    z0_ind=zreq_ind;
    z1_ind=zreq_ind;
    phi0_ind=phireq_ind;
    phi1_ind=phireq_ind;
elseif numel(index_in)==4
    dist.f= h5read(dist_name,'/f',[e_min,p_min,1,1,1],[e_max-e_min+1 p_max-p_min+1 Inf Inf Inf]);
    dist.denf= h5read(dist_name,'/denf');
elseif numel(index_in)==5
    if lz
        zreq_ind=index_in(5);
        zval=dist.z(zreq_ind);
        if ltrim
            dist.f= h5read(dist_name,'/f',[e_min,p_min,1,zreq_ind,1],[e_max-e_min+1 p_max-p_min+1 Inf 1 Inf]);
            dist.denf= h5read(dist_name,'/denf',[1,zreq_ind,1], [Inf 1 Inf]);
            eq=modifyArrays(eq,eq.fields.nphi);
            dist=modifyArrays(dist,dist.nphi);
            phi1_ind=dist.nphi;
        else
            dist.f= h5read(dist_name,'/f',[e_min,p_min,1,zreq_ind,1],[e_max-e_min+1 p_max-p_min+1 Inf 1 Inf]);
            dist.denf= h5read(dist_name,'/denf',[1,zreq_ind,1], [Inf 1 Inf]);
        end
    else
        if dist.nphi>1
            phireq_ind=index_in(5);%
            dist.f= h5read(dist_name,'/f',[e_min,p_min,1,1,phireq_ind],[e_max-e_min+1 p_max-p_min+1 Inf Inf 1]);
            dist.denf= h5read(dist_name,'/denf',[1,1,phireq_ind], [Inf Inf 1]);
        else
            phireq_ind=1;%
            phi1_ind=1;
            dist.f= h5read(dist_name,'/f',[e_min,p_min,1,1],[e_max-e_min+1 p_max-p_min+1 Inf Inf]);
            dist.denf= h5read(dist_name,'/denf',[1,1], [Inf Inf]);
        end    
    end
elseif numel(index_in)==6
    [~,r0_ind]=min(abs(eq.fields.r-index_in(1)));
    [~,z0_ind]=min(abs(eq.fields.z-index_in(5)));
    [~,r1_ind]=min(abs(eq.fields.r-index_in(2)));
    [~,z1_ind]=min(abs(eq.fields.z-index_in(6)));
    if eq.fields.nphi > 1
        [~,phi0_ind]=min(abs(eq.fields.phi-index_in(3)));
        [~,phi1_ind]=min(abs(eq.fields.phi-index_in(4)));
        fprintf('R=%.2f, Phi=%.2f, Z=%.2f\n',eq.fields.r(r0_ind),eq.fields.phi(phi0_ind),eq.fields.z(z0_ind))
        dist.f= h5read(dist_name,'/f',[1,1,r0_ind,z0_ind,phi0_ind],[Inf Inf r1_ind-r0_ind+1 1+z1_ind-z0_ind 1+phi1_ind-phi0_ind]);
        dist.denf= h5read(dist_name,'/denf',[r0_ind,z0_ind,phi0_ind],[r1_ind-r0_ind+1 1+z1_ind-z0_ind 1+phi1_ind-phi0_ind]);
    else
        fprintf('R=%.2f,  Z=%.2f\n',eq.fields.r(r0_ind),eq.fields.z(z0_ind))
        dist.f= h5read(dist_name,'/f',[1,1,r0_ind,z0_ind],[Inf Inf 1 1]);
        dist.denf= h5read(dist_name,'/denf',[r0_ind,z0_ind]);
    end

end
end