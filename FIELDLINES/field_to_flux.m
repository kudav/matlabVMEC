function output = field_to_flux(data,k,mpol,ntor,lasym,varargin)
%FIELD_TO_FLUX(data,k,mpol,ntor,lasym,options) Creates a flux surface from field line data
%   The FIELD_TO_FLUX routine uses FMINSEARCH to fit a set of Fourier
%   harmonics to a set of Poincare data, as generated by the FIELDLINES
%   code.


xm = [];
xn = [];
rmnc = [];
zmns = [];

if ~strcmp(data.datatype,'FIELDLINES')
    output = -1;
    return;
end

% Handle varargin
if nargin > 5
    i=1;
    while i <= length(varargin)
        switch varargin{i}
            case 'xm'
                i=i+1;
                xm_in = varargin{i};
            case 'xn'
                i=i+1;
                xn_in = varargin{i};
            case 'rmnc'
                i=i+1;
                rmnc_in = varargin{i};
            case 'zmns'
                i=i+1;
                zmns_in = varargin{i};
        end
        i = i+1;
    end
end


% Construct Fourier modes
i=1;
for n=0:ntor
    xm(i) = 0;
    xn(i) = -n.*data.nfp;
    i=i+1;
end
for m = 1:mpol
    for n=-ntor:ntor
        xm(i)=m;
        xn(i)=n.*data.nfp;
        i = i + 1;
    end
end
mnmax = length(xm);

% Initialize Rmn and Zmn
clf;
if isempty(rmnc_in)
    rmnc = zeros(mnmax,1);
    zmns = zeros(mnmax,1);
    dex = and((xm == 0),(xn == 0));
    rmnc(dex,1) = mean(data.R_lines(1,:));
    zmns(dex,1) = 0.0;
    dex = and((xm == 1),(xn == 0));
    drho = mean(abs(data.R_lines(k,1:data.npoinc:end)-data.R_lines(1,1:data.npoinc:end)));
    dzed = mean(abs(data.Z_lines(k,1:data.npoinc:end)-data.Z_lines(1,1:data.npoinc:end)));
    rmnc(dex,1) = drho;
    zmns(dex,1) = dzed;
    drho2=1.5;
    for i=1:mnmax
        if (rmnc(i,1) ~= 0), continue; end
        if (zmns(i,1) ~= 0), continue; end
        rmnc(i,1) = (rand-0.5)*drho/(2+drho2^xm(i)+drho2^abs(xn(i)./data.nfp));
        zmns(i,1) = (rand-0.5)*dzed/(2+drho2^xm(i)+drho2^abs(xn(i)./data.nfp));
    end
else
    rmnc = zeros(mnmax,1);
    zmns = zeros(mnmax,1);
    for i=1:mnmax
        dex = and(xm_in == xm(i), xn_in == xn(i));
        if isempty(dex), continue; end
        temp = rmnc_in(dex);
        if ~isempty(temp)
            rmnc(i,1) = temp;
            temp = zmns_in(dex);
            zmns(i,1) = temp;
        end
    end
end

% Initialize
dex = data.R_lines(k,1:720) > 0;
xtarg = data.X_lines(k,dex);
ytarg = data.Y_lines(k,dex);
ztarg = data.Z_lines(k,dex);
ptarg = data.PHI_lines(k,dex);
i = find(ptarg > 2*pi,1,'first');
zeta = [ptarg(1:i-1) 2*pi];
npoinc = data.npoinc;
mn = mnmax;
nt = 2.*length(xtarg)/(i);
nt = 128;
%nt=round(data.nsteps./data.npoinc);
%nt = nt * 2;
theta =0:2*pi/double(nt-1):2*pi;
%zeta  =0:2*pi/double(data.npoinc):2*pi;

% Define the problem
%x0=[rmnc.^(1+xm');zmns(2:end).^(1+xm(2:end)')];
x0=[rmnc;zmns(2:end)];
opts=optimset('MaxIter',10000,...
    'TolFun',1.0E-8,...
    'Display','iter',...
    'OutputFcn',@(x,optimvalues,states) pltfun(x,optimvalues,states,mn,xm,xn,theta,zeta,xtarg,ytarg,ztarg,npoinc));

% Run the problem (fminsearch)
x=fminsearch(@(x) minfun(x,mn,xm,xn,theta,zeta,xtarg,ytarg,ztarg),x0,opts);

% Run the problem (Levenberg-Marquardt)
% t0 = zeros(1,length(xtarg));
% y0 = zeros(1,length(xtarg));
% weight = ones(1,length(xtarg));
% dp = 0.001;
% p_min = -10.*ones(1,length(x0));
% p_max = 10.*ones(1,length(x0));
% c = [ mn length(theta) length(zeta) length(xtarg) xm xn theta zeta xtarg ytarg ztarg];
% 
% [p, redX2, sigma_p,sigma_y,corr_p,R_sq,cvg_hst] = ...
%     lm(@(t,p,c) lm_func(t,p,c) ,x0,t0,y0,weight,dp,p_min,p_max,c);
% 
% x=p;
% Return the output
output.rmnc=x(1:mnmax);
output.zmns(1) = 0;
output.zmns(2:mnmax)=x(mnmax+1:end);
output.xm = xm;
output.xn = xn;
return
end

function F = minfun(x,mn,xm,xn,theta,zeta,xtarg,ytarg,ztarg)
%rmn = x(1:mn).^1./(1+xm');
%zmn(2:mn,1) = x(mn+1:end).^1./(1+xm(2:end)');
rmn = x(1:mn);
zmn(2:mn,1) = x(mn+1:end);
zmn(1,1) = 0;
r   = squeeze(cfunct(theta,zeta,rmn,xm,xn));
z2   = squeeze(sfunct(theta,zeta,zmn,xm,xn));
for i=1:size(r,2)
    x2(:,i) = r(:,i).*cos(zeta(i));
    y2(:,i) = r(:,i).*sin(zeta(i));
end
F = 0;
for i=1:length(xtarg)
    dx = x2-xtarg(i);
    dy = y2-ytarg(i);
    dz = z2-ztarg(i);
    dl = dx.*dx+dy.*dy+dz.*dz;
    F = F+sqrt(min(min(dl)));
end
F = F./length(xtarg);
return
end

function F = lm_func(x,p,c)
%rmn = x(1:mn).^1./(1+xm');
%zmn(2:mn,1) = x(mn+1:end).^1./(1+xm(2:end)');
mn = c(1);
ntheta = c(2);
nzeta = c(3);
ntarg = c(4);
istart = 5; iend = 5+mn-1; xm = c(istart:iend);
istart = iend+1; iend = istart+mn-1; xn = c(istart:iend);
istart = iend+1; iend = istart+ntheta-1; theta = c(istart:iend);
istart = iend+1; iend = istart+nzeta-1; zeta = c(istart:iend);
istart = iend+1; iend = istart+ntarg-1; xtarg = c(istart:iend);
istart = iend+1; iend = istart+ntarg-1; ytarg = c(istart:iend);
istart = iend+1; iend = istart+ntarg-1; ztarg = c(istart:iend);
rmn = p(1:mn);
zmn(2:mn,1) = p(mn+1:end);
zmn(1,1) = 0;
r   = squeeze(cfunct(theta,zeta,rmn,xm,xn));
z2   = squeeze(sfunct(theta,zeta,zmn,xm,xn));
for i=1:size(r,2)
    x2(:,i) = r(:,i).*cos(zeta(i));
    y2(:,i) = r(:,i).*sin(zeta(i));
end
F = zeros(1,length(xtarg));
for i=1:length(xtarg)
    dx = x2-xtarg(i);
    dy = y2-ytarg(i);
    dz = z2-ztarg(i);
    dl = dx.*dx+dy.*dy+dz.*dz;
    F(i) = sqrt(min(min(dl)));
end
return
end


function stop = pltfun(x,optimvalues,states,mn,xm,xn,theta,zeta,xtarg,ytarg,ztarg,npoinc)
stop = false;
%rmn = x(1:mn).^1./(1+xm');
%zmn(2:mn,1) = x(mn+1:end).^1./(1+xm(2:end)');
rmn = x(1:mn);
zmn(2:mn,1) = x(mn+1:end);
zmn(1,1) = 0;
r   = cfunct(theta,zeta,rmn,xm,xn);
z2   = sfunct(theta,zeta,zmn,xm,xn);
subplot(2,2,1);
cla;
plot(r(1,:,1),z2(1,:,1),'r','Linewidth',2.0);
hold on;
rtarg = sqrt(xtarg.^2+ytarg.^2);
plot(rtarg(1:npoinc:end),ztarg(1:npoinc:end),'.k');
hold off;
axis equal;
subplot(2,2,3);
cla;
nfp = min(xn(xn>0));
nz2=round(length(zeta)./(2*nfp))+1;
plot(r(1,:,nz2),z2(1,:,nz2),'r','Linewidth',2.0);
hold on;
nz2=npoinc/2+1;
rtarg = sqrt(xtarg.^2+ytarg.^2);
plot(rtarg(nz2:npoinc:end),ztarg(nz2:npoinc:end),'.k');
hold off;
axis equal;
subplot(2,2,2);
cla;
bar(x);
hold on;
plot([1 1].*mn,ylim,'r');
subplot(2,2,4);
cla;
isotoro(r,z2,zeta,1);
hold on;
plot3(xtarg,ytarg,ztarg,'.k');
axis equal;
hold off;
drawnow;
return
end



