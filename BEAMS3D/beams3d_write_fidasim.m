function [f, denf] = beams3d_write_fidasim(data, name,varargin)
% BEAMS3D_WRITE_FIDASIM produces the FIDASIM input files after a run of BEAMS3D.
%
% For new versions of BEAMS3D, the same functionality can be achieved with the 
% -fidasim flag. The distribution function and all quantities are output in the 
% standard BEAMS3D cylindrical grid as standard.
%
% Optional parameters can be used to customize the grid spacing or set the axis 
% precisely:
%   - 'n', nr, nphi, nz, nE, np: Change the grid spacing.
%   - 'axis', raxis, phiaxis, zaxis: Set the axis precisely.
%
% The energy range is set to 0-Emax from the maximum particle velocity. Flow 
% velocities are set to 0, and the electric field is calculated from the gradient 
% of POT_ARR. Optionally, the function outputs the fast ion distribution and density 
% as variables.
%
% Example usage:
%   data = read_beams3d('test.h5');
%   beams3d_write_fidasim(data, 'fidasim_test', 'n', nr, nphi, nz, nE, np);
%   beams3d_write_fidasim(data, 'fidasim_test', 'axis', raxis, phiaxis, zaxis);
%
% Inputs:
%   - data: The BEAMS3D data structure obtained from read_beams3d.
%   - output_filename: The name of the output file for FIDASIM input.
%
% Optional Parameters:
%   - 'n': Change the grid spacing.
%       nr: Radial grid points.
%       nphi: Azimuthal grid points.
%       nz: Axial grid points.
%       nE: Energy grid points.
%       np: Pitch angle grid points.
%   - 'axis': Set the axis precisely.
%       raxis: Radial axis values.
%       phiaxis: Azimuthal axis values.
%       zaxis: Axial axis values.



ec  = 1.60217662E-19; % electron charge [C]
amu = 1.66053906660E-27; % Dalton [kg]
lmovie = 0;
lrecalc=0;
inputs=0;
beam_dex = 1:data.nbeams;

filename_dist = [name,'_distribution.h5'];
filename_eq = [name,'_equilibrium.h5'];
time = 0.0;
mass = data.mass(1);

if isfile(filename_dist)
    disp('Found old distribution file. Overwriting...')
    delete(filename_dist);
end
if isfile(filename_eq)
    disp('Found old equilibrium file. Overwriting...')
    delete(filename_eq);
end

%Get grid points
% Return an RPZ shaped array;
raxis   = data.raxis';%linspace(min(data.raxis), max(data.raxis), 70);
zaxis   = data.zaxis';%linspace(min(data.zaxis), max(data.zaxis), 71);
paxis   = data.phiaxis';%linspace(0, 2*pi, 15);
nE = data.ns_prof4;
np = data.ns_prof5;
dr=raxis(2)-raxis(1);
dz=zaxis(2)-zaxis(1);
dphi=paxis(2)-paxis(1);
% Handle varargin
if ~isempty(varargin)
    i=1;
    while i <= length(varargin)
        switch varargin{i}
            case{'axis'}
                i = i + 1;
                raxis = varargin{i};
                if size(raxis,2)==1
                    raxis = raxis';
                end
                i = i + 1;
                paxis = varargin{i};
                if size(paxis,2)==1
                    paxis = paxis';
                end
                i = i + 1;
                zaxis = varargin{i};
                if size(zaxis,2)==1
                    zaxis = zaxis';
                end
            case{'n'}
                i = i + 1;
                dr=raxis(2)-raxis(1);
                raxis   =  ((1:varargin{i})-1).*(raxis(end)-raxis(1)+dr)/varargin{i} + raxis(1)-dr/2;%!Lower grid edges%linspace(min(raxis), max(raxis), varargin{i});
                i = i + 1;
                dphi=paxis(2)-paxis(1);
                paxis   = ((1:varargin{i})-1).*(paxis(end)-paxis(1)+dphi)/varargin{i}  + paxis(1)-dphi/2;%linspace(min(data.phiaxis), max(data.phiaxis), varargin{i});
                i = i + 1;
                dz=zaxis(2)-zaxis(1);
                zaxis   = ((1:varargin{i})-1).*(zaxis(end)-zaxis(1)+dz)/varargin{i} + zaxis(1)-dz/2;
                %zaxis   = linspace(min(data.zaxis), max(data.zaxis), varargin{i});
                i = i + 1;
                nE = varargin{i};
                i = i + 1;
                np = varargin{i};
            case{'movie'}
                lmovie = 1;
            case{'recalc_dist'}
                lrecalc=1;
                i=i+1;
                inputs=varargin{i};
                i=i+1;
                type=varargin{i};
            case{'beams'}
                i = i+1;
                beam_dex = varargin{i};
            otherwise
                disp(['Unrecognized Option: ' varargin{i}]);
                return
        end
        i = i + 1;
    end
end



Emax = ceil(0.5.*mass.*data.partvmax.^2./ec);
%Eaxis   = 0.5*max(Emax)/nE:max(Emax)/nE:max(Emax);%linspace(0,Emax,nE);%%0:10E3:100E3;
%pitchaxis = -1+1/np:2/np:1;%linspace(-1,1,np);%-1:0.1:1;
%Same axis as from BEAMS3D interface:
Eaxis = ((1:nE) - 0.5) / nE * Emax;
pitchaxis = ((1:np) - 0.5) / np * 2 - 1;

[R,P,Z,E,PITCH] = ndgrid(raxis,paxis,zaxis,Eaxis,pitchaxis);

[R3,P3,Z3] = ndgrid(raxis,paxis,zaxis);
[Rd,Pd,Zd] = ndgrid(data.raxis,data.phiaxis,data.zaxis);
P3_mod=mod(P3,max(data.phiaxis));
%nsave = size(R);
ntotal = numel(R);
R = reshape(R,[1 ntotal]);
P = reshape(P,[1 ntotal]);
Z = reshape(Z,[1 ntotal]);
E = reshape(E,[1 ntotal]);
PITCH = reshape(PITCH,[1 ntotal]);

%Convert to cm
raxis   = raxis.*100;
%R3 = R3 .*100;
zaxis   = zaxis.*100;
%Z3 = Z3 .*100;
Eaxis = Eaxis./1000;
%E = E/1000;

if lrecalc
    if numel(inputs)==0
        inputs{1}=3; %Indices for orbit source (S(rho,xi)), for npoinc=5 combined runs, this should be 3 or 3:4
        inputs{2} = 75; %npitch
        inputs{3} = 50; %nv
        inputs{4}=50; %norder, Order up to which to calculate the legendre polynomials
        inputs{5}=[1.0 1.0];%mi/ Ai, in amu
        inputs{6}=[1.0 1.0];%Zi in ec
    end
    [data, rho, pitchaxis, Eaxis, f2D,~,~,~] = mRabbit(data, inputs{1}, inputs{2}, inputs{3}, inputs{4},inputs{5},inputs{6},type);
    f2D=squeeze(sum(f2D,3,'omitnan'));

    n=squeeze(trapz(pitchaxis,trapz(Eaxis,f2D,1),2));
    F = griddedInterpolant(Rd-dr/2,Pd-dphi/2,Zd-dz/2,sqrt(data.S_ARR), 'spline');
    RHO_ARR = F(R3, P3_mod, Z3);

    %RHO_ARR=sqrt(data.S_ARR(:,1:5,:));

    denf=interp1(rho,n,RHO_ARR(:),'linear',0);
    denf=reshape(denf,size(RHO_ARR));
    denf=permute(denf,[1 3 2]);

    f=interp1(rho,reshape(f2D,[inputs{3}*inputs{2} numel(rho)])',RHO_ARR(:),'linear',0);
    f=reshape(f,[size(RHO_ARR), inputs{3} inputs{2}]);
    f=permute(f,[4 5 1 3 2]);
    %Quick fix
    f(isnan(f))=0;
    f(isinf(f))=0;
else
    f=beams3d_getdistrpzEpitch(data,R,P,Z,E,PITCH);
    f = f./(1000*1E6); %ev -> keV, m^-3 -> cm^-3
    f = sum(f(beam_dex,:),1); % Sum over beamlines
    f = reshape(f,[numel(raxis), numel(paxis), numel(zaxis), numel(Eaxis), numel(pitchaxis)]);
    f= permute(f,[4, 5, 1, 3, 2]); %Align with FIDASIM axis order

    %Quick fix
    %f(isnan(f))=0;
    %f(isinf(f))=0;

    denf = squeeze(trapz(Eaxis, f,1));
    denf = squeeze(trapz(pitchaxis, denf,1));%Integration in velocity space
end


F = griddedInterpolant(Rd,Pd,Zd, data.B_R, 'spline');
br = F(R3, P3_mod, Z3);
F = griddedInterpolant(Rd,Pd,Zd, data.B_PHI, 'spline');
bt = F(R3, P3_mod, Z3);
F = griddedInterpolant(Rd,Pd,Zd, data.B_Z, 'spline');
bz = F(R3, P3_mod, Z3);

vr = zeros(size(br));
if isfield(data,'VTOR_ARR')
    F = griddedInterpolant(Rd,Pd,Zd, data.VTOR_ARR*100, 'spline');
    vt = F(R3, P3_mod, Z3);
else
    vt = vr;
end
vz = vr;

[er, et, ez] = gradient(data.POT_ARR,mean(diff(raxis)),mean(diff(paxis)),mean(diff(zaxis))); %first output corresponds to gradient along 2nd dimension???
F = griddedInterpolant(Rd,Pd,Zd, er, 'spline');
er = F(R3, P3_mod, Z3);
F = griddedInterpolant(Rd,Pd,Zd, et, 'spline');
et = F(R3, P3_mod, Z3);
F = griddedInterpolant(Rd,Pd,Zd, ez, 'spline');
ez = F(R3, P3_mod, Z3);

F = griddedInterpolant(Rd,Pd,Zd, data.NE/1e6, 'spline');
ne = F(R3, P3_mod, Z3);
F = griddedInterpolant(Rd,Pd,Zd, data.TE/1000, 'spline');
te = F(R3, P3_mod, Z3);
F = griddedInterpolant(Rd,Pd,Zd, data.TI/1000, 'spline');
ti = F(R3, P3_mod, Z3);
F = griddedInterpolant(Rd,Pd,Zd, data.ZEFF_ARR, 'spline');
zeff = F(R3, P3_mod, Z3);
%denn = zeros(size(zeff));


write_new_data_with_att(filename_dist,'/type',1,'description', 'Distribution type: 1="Guiding Center Density Function", 2="Guiding Center Monte Carlo", 3="Full Orbit Monte Carlo"')
write_new_data_with_att(filename_dist,'/f',f,{'description', 'units'}, {'Distribution Function (nenergy_fida,npitch_fida,nr_fida,nz_fida,nphi_fida)', 'part/(cm^3 keV)'})
write_new_data_with_att(filename_dist,'/denf',denf,{'description', 'units'}, {'Fast-ion density (nr_fida,nz_fida,nphi_fida)','part/(cm^3)'})

write_new_data_with_att(filename_dist,'/nenergy',numel(Eaxis),'description','Number of energy values');
write_new_data_with_att(filename_dist,'/npitch',numel(pitchaxis), 'description','Number of pitch values');
write_new_data_with_att(filename_dist,'/nr',numel(raxis),'description','Number of R values');
write_new_data_with_att(filename_dist,'/nphi',numel(paxis),'description','Number of Phi values');
write_new_data_with_att(filename_dist,'/nz', numel(zaxis),'description','Number of Z values');
write_new_data_with_att(filename_dist,'/time',time,{'description', 'units'},{'Distribution time', 's'});
write_new_data_with_att(filename_dist,'/r',raxis,{'description', 'units'},{'Radius','cm'});
write_new_data_with_att(filename_dist,'/phi',paxis,{'description', 'units'},{'Toroidal angle','rad'});
write_new_data_with_att(filename_dist,'/z',zaxis,{'description', 'units'},{'Z','cm'});
write_new_data_with_att(filename_dist,'/energy',Eaxis,{'description', 'units'},{'Energy array','keV'});
write_new_data_with_att(filename_dist,'/pitch',pitchaxis,{'description', 'units'},{'Pitch array','-'});
write_new_data_with_att(filename_dist,'/r2d',repmat(raxis',1,numel(zaxis)),{'description', 'units'},{'R grid: R(r,z)','cm'});
write_new_data_with_att(filename_dist,'/z2d',repmat(zaxis,numel(raxis),1),{'description', 'units'},{'Z grid: Z(r,z)','cm'});


% EQUILIBRIUM FILE
%Fields
write_new_data_with_att(filename_eq,'/fields/nr',numel(raxis),'description','Number of R values');
write_new_data_with_att(filename_eq,'/fields/nphi',numel(paxis),'description','Number of Phi values');
write_new_data_with_att(filename_eq,'/fields/nz', numel(zaxis),'description','Number of Z values');
write_new_data_with_att(filename_eq,'/fields/time',time,'description','Distribution time');

write_new_data_with_att(filename_eq,'/fields/r',raxis,{'description', 'units'},{'Radius','cm'});
write_new_data_with_att(filename_eq,'/fields/phi',paxis,{'description', 'units'},{'Toroidal angle','rad'});
write_new_data_with_att(filename_eq,'/fields/z',zaxis,{'description', 'units'},{'Z','cm'});

write_new_data_with_att(filename_eq,'/fields/r2d',repmat(raxis',1,numel(zaxis)),{'description', 'units'},{'R grid: R(r,z)','cm'});
write_new_data_with_att(filename_eq,'/fields/z2d',repmat(zaxis,numel(raxis),1),{'description', 'units'},{'Z grid: Z(r,z)','cm'});

write_new_data_with_att(filename_eq,'/fields/mask',int32(ones(numel(raxis),numel(zaxis), numel(paxis))),'description','Boolean mask that indicates where the fields are well defined');



write_new_data_with_att(filename_eq,'/fields/br',permute(br,[1,3,2]),{'description', 'units'},{'Magnetic field in the r-direction: Br(r,z,phi)','T'});
write_new_data_with_att(filename_eq,'/fields/bt',permute(bt,[1,3,2]),{'description', 'units'},{'Magnetic field in the toroidal phi-direction: Bt(r,z,phi)','T'});
write_new_data_with_att(filename_eq,'/fields/bz',permute(bz,[1,3,2]),{'description', 'units'},{'Magnetic field in the z-direction: Bz(r,z,phi)','T'});


write_new_data_with_att(filename_eq,'/fields/er',permute(-er,[1,3,2]),{'description', 'units'},{'Electric field in the r-direction: Er(r,z,phi)','V/m'});
write_new_data_with_att(filename_eq,'/fields/et',permute(-et,[1,3,2]),{'description', 'units'},{'Electric field in the toroidal phi-direction: Et(r,z,phi)','V/m'});
write_new_data_with_att(filename_eq,'/fields/ez',permute(-ez,[1,3,2]),{'description', 'units'},{'Electric field in the toroidal phi-direction: Ez(r,z,phi)','V/m'});

%Plasma
write_new_data_with_att(filename_eq,'/plasma/nr',numel(raxis),'description','Number of R values');
write_new_data_with_att(filename_eq,'/plasma/nphi',numel(paxis),'description','Number of Phi values');
write_new_data_with_att(filename_eq,'/plasma/nz', numel(zaxis),'description','Number of Z values');
write_new_data_with_att(filename_eq,'/plasma/time',time,'description','Distribution time');

write_new_data_with_att(filename_eq,'/plasma/r',raxis,{'description', 'units'},{'Radius','cm'});
write_new_data_with_att(filename_eq,'/plasma/phi',paxis,{'description', 'units'},{'Toroidal angle','rad'});
write_new_data_with_att(filename_eq,'/plasma/z',zaxis,{'description', 'units'},{'Z','cm'});

write_new_data_with_att(filename_eq,'/plasma/r2d',repmat(raxis',1,numel(zaxis)),{'description', 'units'},{'R grid: R(r,z)','cm'});
write_new_data_with_att(filename_eq,'/plasma/z2d',repmat(zaxis,numel(raxis),1),{'description', 'units'},{'Z grid: Z(r,z)','cm'});

write_new_data_with_att(filename_eq,'/plasma/mask',int32(ones(numel(raxis),numel(zaxis), numel(paxis))),'description','Boolean mask that indicates where the fields are well defined');


% write_new_data_with_att(filename_eq,'/plasma/er',zeros(numel(raxis),numel(zaxis)),{'description', 'units'},{'Bulk plasma flow in the r-direction: Vr(r,z,phi)','cm/s'});
% write_new_data_with_att(filename_eq,'/plasma/et',zeros(numel(raxis),numel(zaxis)),{'description', 'units'},{'Bulk plasma flow in the toroidal phi-direction: Vt(r,z,phi)','cm/s'});
% write_new_data_with_att(filename_eq,'/plasma/ez',zeros(numel(raxis),numel(zaxis)),{'description', 'units'},{'Bulk plasma flow in the toroidal phi-direction: Vz(r,z,phi)','cm/s'});
write_new_data_with_att(filename_eq,'/plasma/vr',permute(vr,[1,3,2]),{'description', 'units'},{'Bulk plasma flow in the r-direction: Vr(r,z,phi)','cm/s'});
write_new_data_with_att(filename_eq,'/plasma/vt',permute(vt,[1,3,2]),{'description', 'units'},{'Bulk plasma flow in the toroidal phi-direction: Vt(r,z,phi)','cm/s'});
write_new_data_with_att(filename_eq,'/plasma/vz',permute(vz,[1,3,2]),{'description', 'units'},{'Bulk plasma flow in the toroidal phi-direction: Vz(r,z,phi)','cm/s'});



write_new_data_with_att(filename_eq,'/plasma/dene',permute(ne,[1,3,2]),{'description', 'units'},{'Electron Number Density: Dene(r,z, phi)','cm^-3'});
write_new_data_with_att(filename_eq,'/plasma/te',permute(te,[1,3,2]),{'description', 'units'},{'Electron Temperature: Te(r,z,phi)','keV'});
write_new_data_with_att(filename_eq,'/plasma/ti',permute(ti,[1,3,2]),{'description', 'units'},{'Ion Temperature: Ti(r,z,phi)','keV'});
write_new_data_with_att(filename_eq,'/plasma/zeff',permute(zeff,[1,3,2]),{'description', 'units'},{'Effective Nuclear Charge: Zeff(r,z,phi)','-'});
%write_new_data_with_att(filename_eq,'/plasma/denn',permute(denn,[1,3,2]),{'description', 'units'},{'Neutral density: denn(r,z,phi)','cm^-3'});

if lmovie
    v= VideoWriter(['output_' name '.avi']);
    open(v);
    f = figure;
    for i =1:numel(paxis)
        subplot(3,3,1)
        pixplot(raxis,zaxis,squeeze(br(:,i,:)))
        title('BR')
        xlabel('R [m]')
        ylabel('Z [m]')
        subplot(3,3,2)
        pixplot(raxis,zaxis,squeeze(bt(:,i,:)))
        title('BT')
        xlabel('R [m]')
        ylabel('Z [m]')
        subplot(3,3,3)
        pixplot(raxis,zaxis,squeeze(bz(:,i,:)))
        title('BZ')
        xlabel('R [m]')
        ylabel('Z [m]')
        subplot(3,3,4)
        pixplot(raxis,zaxis,squeeze(er(:,i,:)))
        title('ER')
        xlabel('R [m]')
        ylabel('Z [m]')
        subplot(3,3,5)
        pixplot(raxis,zaxis,squeeze(et(:,i,:)))
        title('ET')
        xlabel('R [m]')
        ylabel('Z [m]')
        subplot(3,3,6)
        pixplot(raxis,zaxis,squeeze(ez(:,i,:)))
        title('EZ')
        xlabel('R [m]')
        ylabel('Z [m]')
        subplot(3,3,7)
        pixplot(raxis,zaxis,squeeze(ne(:,i,:)))
        title('NE')
        xlabel('R [m]')
        ylabel('Z [m]')
        subplot(3,3,8)
        pixplot(raxis,zaxis,squeeze(te(:,i,:)))
        title('TE')
        xlabel('R [m]')
        ylabel('Z [m]')
        subplot(3,3,9)
        pixplot(raxis,zaxis,squeeze(ti(:,i,:)))
        title('TI')
        xlabel('R [m]')
        ylabel('Z [m]')

        frame = getframe(f);
        writeVideo(v,frame);
    end
end

end

function write_new_data_with_att(filename,loc,data,attrs,vals)
if sum(size(data)==1)==1
    h5create(filename,loc, length(data)); %Initialize only one dimension for vectors!
else
    h5create(filename,loc, size(data));
end
h5write(filename,loc,squeeze(data));
if iscell(attrs)
    for i = 1:numel(attrs)
        h5writeatt(filename,loc,attrs{i},vals{i})
    end
else
    h5writeatt(filename,loc,attrs,vals)
end
end