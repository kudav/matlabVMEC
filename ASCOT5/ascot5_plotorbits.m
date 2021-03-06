function fig = ascot5_plotorbits(a5file,runid,varargin)
%ASCOT5_PLOTORBITS plots the orbit data from an ASCOT5 file
%   The ASCOT5_PLOTORBITS routine plots the 'orbits' field of an ASCOT5
%   HDF5 file. It takes the ASCOT5 filename, and runid as input. If an
%   empty array is passed then the active id for orbits is used. Optional
%   arguments control what type of plot is generated.
%
%   Options arguments:
%       'parts':    Select subset of Markers by ID
%       'xyz':      XYZ orbit plot (default)
%       'flux':     Rho-theta polar plot
%       'fluxperp': Rho-theta xy-plot
%       'pll':      Parallel-Perpendicular Momentum plot
%       'rz':       Plot in R-Z Plane
%       'movie':    Make a movie of the plot
%
%   Example:
%       a5file='ascot5_test.h5';
%       orbitid=0396210459;
%       fig=ascot5_plotorbits(a5file,orbitid);
%       fig=ascot5_plotorbits(a5file,[]); %Use active ID
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

% Defaults
amu = 1.66053906660E-27;
plottype='xyz';
parts = [];
lmovie = 0;
pspace = 0;
filename = 'default';
startframe = 100; %Number of Start-Frame - ntail
endframe = 20000; %Number of last frame
skipframes = 30; %How many data points to skip between images
ntail = 20;

% Handle varargin
if nargin > 2
    i=1;
    while i <= nargin-2
        switch varargin{i}
            case{'xyz','flux','vll','pll', 'rz', 'fluxperp'}
                plottype=varargin{i};
            case 'parts'
                i = i+1;
                parts=varargin{i};
            case 'movie'
                i = i + 1;
                filename = varargin{i};
                lmovie = 1;
            case 'ntail'
                i = i+1;
                ntail = varargin{i};
            case {'pspace', 'vspace'}
                pspace = 1;
            otherwise
                disp(['Unrecognized Option: ' varargin{i}]);
                return
        end
        i = i + 1;
    end
end

fig=[];
% Check for file
if ~isfile(a5file)
    disp(['ERROR: ' a5file ' file not found!']);
    return;
end

% Use active run
if isempty(runid)
    try
        runid=h5readatt(a5file,'/results','active');
        disp(['  Using runid: ' runid]);
    catch
        runid=[];
        return;
    end
end

% Set run path
path = ['/results/run_' num2str(runid,'%10.10i') '/orbit'];

% Get IDS and time for sorting
ids = double(h5read(a5file,[path '/ids']));
time = double(h5read(a5file,[path '/mileage']));
if isempty(parts)
    parts = unique(ids);
end


% get array sizes

nsteps = round(length(ids)./max(ids));
endframe = nsteps;
dex = ismember(ids,parts);
ids = ids(dex);
time = time(dex);

npart = length(unique(parts));

time = reshape(time,[nsteps npart]);

% create particle sorting array
[time,idex] = sort(time);


%Calculate Orbits
% Extract position
rho = h5read(a5file,[path '/rho']);
theta = h5read(a5file,[path '/theta']);
phi = h5read(a5file,[path '/phi']);

r = h5read(a5file,[path '/r']);
z = h5read(a5file,[path '/z']);
x = r.*cos(phi);
y = r.*sin(phi);

pr = h5read(a5file,[path '/pr']);
pz = h5read(a5file,[path '/pz']);
pphi = h5read(a5file,[path '/pphi']);

br = h5read(a5file,[path '/br']);
bz = h5read(a5file,[path '/bz']);
bphi = h5read(a5file,[path '/bphi']);

% get pll
b=sqrt(br.*br+bphi.*bphi+bz.*bz);
pll = (pr.*br+pphi.*bphi+pz.*bz)./b;
p2  = pr.*pr+pphi.*pphi+pz.*pz;
pperp = sqrt(p2-pll.*pll);

b = b(dex);
r = r(dex);

x = x(dex);
y = y(dex);
z = z(dex);
rho = rho(dex);
theta = theta(dex);
phi = phi(dex);
pll = pll(dex);
pperp = pperp(dex);


% Reshape by particle index
b = reshape(b,[nsteps npart]);
r = reshape(r,[nsteps npart]);
x = reshape(x,[nsteps npart]);
y = reshape(y,[nsteps npart]);
z = reshape(z,[nsteps npart]);
rho = reshape(rho,[nsteps npart]);
theta = reshape(theta,[nsteps npart]);
phi = reshape(phi,[nsteps npart]);
pll = reshape(pll,[nsteps npart]);
pperp = reshape(pperp,[nsteps npart]);

% Reorder in time
for i=1:npart
    b(:,i) = b(idex(:,i),i);
    r(:,i) = r(idex(:,i),i);
    x(:,i) = x(idex(:,i),i);
    y(:,i) = y(idex(:,i),i);
    z(:,i) = z(idex(:,i),i);
    rho(:,i) = rho(idex(:,i),i);
    theta(:,i) = theta(idex(:,i),i);
    phi(:,i) = phi(idex(:,i),i);
    pll(:,i) = pll(idex(:,i),i);
    pperp(:,i) = pperp(idex(:,i),i);
end
%----------Calculate Orbits

%Prepare Figure
fig = figure('Position',get(0,'ScreenSize'));%figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
ax1 = gca;
if pspace
    delete(gca);
    ax1 = axes('Position', [0.07 0.18 0.55 0.75]);
    ax2 = axes('Position', [.75 0.18 0.2*2/3.0 0.2]);
    ax3 = axes('Position', [.68 0.58 0.3 0.4]);
end
%linecolors = num2cell(lines(size(x,2)),2);
linecolors = lines(size(x,2));
v = VideoWriter([filename '.mp4'], 'MPEG-4');
open(v);

switch plottype
    case 'xyz'
        modphi = mod(phi,360);
        edges = linspace(0,360,360);
        discphi = discretize(modphi,edges);
        colors = parula(360);
        % Make plot
        %plot3(ax1,x,y,z);
        hold on;
        if lmovie
            for i=startframe+ntail:skipframes:endframe
                
                %s1 = scatter3(ax1,x(i-ntail:i,:),y(i-ntail:i,:),z(i-ntail:i,:));
                s2 = plot3(ax1,x(i-ntail:i,:),y(i-ntail:i,:),z(i-ntail:i,:),'Color', colors(discphi(i),:));
                campos([95.6559   -2.0102   32.5105]);
                camtarget([-1.0351   -0.2400    0.1142]);
                camup([0 0 1]);
                set(gca,'FontSize',24);
                xlabel('X [m]');
                ylabel('Y [m]');
                zlabel('Z [m]');
                axis equal
                caxis(ax1,[0 360]);
                zlim(ax1,[-1.25 1.25]);
                ylim(ax1,[-6 6]);
                xlim(ax1,[-6 6]);
                
                title('Particle Orbits');
                
                frame = getframe(fig); %for writing
                writeVideo(v,frame);
                %delete(s1);
                delete(s2);
            end
        else
            scatter3(ax1,x,y,z,0.1,time);
            set(gca,'FontSize',24);
            xlabel('X [m]');
            ylabel('Y [m]');
            zlabel('Z [m]');
            title('Particle Orbits');
            axis equal;
        end
    case 'flux'
        
        % Make plot
        delete(ax1);
        ax1 = polaraxes;
        modphi = mod(phi,180);
        edges = linspace(0,180,360);
        discphi = discretize(modphi,edges);
        colors = parula(360);
        
        if lmovie
            for i=startframe+ntail:skipframes:endframe
                polarplot(ax1,deg2rad(theta(i-ntail:i,:)),rho(i-ntail:i,:), 'Color', colors(discphi(i),:));
                ax1.RLim = [0 1.5];
                
                %scatter(ax1, deg2rad(theta(i,:),rho(i,:));
                set(gca,'FontSize',24);
                title('Particle Orbits');
                frame = getframe(fig); %for writing
                writeVideo(v,frame);
            end
        else
            polarplot(ax1,deg2rad(theta),rho);
            set(gca,'FontSize',24);
            title('Particle Orbits');
        end
    case 'fluxperp'
        
        % Make plot 
        %modphi = mod(phi,180);
        %edges = linspace(0,180,360);
        %discquant= discretize(modphi,edges);
        %colors = parula(360);
        edges = linspace(min(b, [], 'all'), max(b, [], 'all'), 360);
        
        colors = parula(360);
        discquant = discretize(b,edges);
        poloidalextent = 360.0;
        toroidalextent = 72.0;
        polrevolutions = floor(theta./ poloidalextent); %How many full revolutions have the particles completed
        torrevolutions = floor(phi ./ toroidalextent);
        %phi = phi - torrevolutions(1,:).*toroidalextent; %Reset rotations
        %theta = theta - polrevolutions(1,:).*poloidalextent;
        
        if lmovie
            tail = ntail;
        else
            tail = endframe - startframe;
            startframe = endframe-ntail;
        end
        
        hold(ax1,'on');
        hold(ax2,'on');
        hold(ax3,'on');
        for i=startframe+ntail:skipframes:endframe
            plottheta = theta(i-tail:i,:);
            plotrho = rho(i-tail:i,:);
            plotcolors = colors(discquant(i-tail:i,:),:);
            plot_polrevolutions = [polrevolutions(i-tail,:)' polrevolutions(i,:)']; %How many times does the orbit wrap around in current plot for each particle
            plot_torrevolutions = [torrevolutions(i-tail,:)' torrevolutions(i,:)'];
            %                 h1 = plot(ax1,plotx,ploty);
            %                 set(h1, {'color'}, linecolors);
            %                 hold(ax1,'on');
            
            %if sum(plotrevolutions) %replot when part of orbit goes offscreen
            for j = 1:numel(parts) %Plot all frames where particle goes over plot area boundary
                if plot_polrevolutions(j, 1)>plot_polrevolutions(j, 2)
                    plot_polrevolutions = fliplr(plot_polrevolutions);
                end
                if plot_torrevolutions(j,1)>plot_torrevolutions(j,2)
                    plot_torrevolutions = fliplr(plot_torrevolutions);
                end
                for k = plot_polrevolutions(j, 1):plot_polrevolutions(j, 2)%Frames for poloidal crossings
                    plotx = plottheta(:,j)-poloidalextent*k;
                    ploty = plotrho(:,j);
                    scattercolors = plotcolors((j-1)*size(plottheta,1)+1:(j*size(plottheta,1)),:);
                    dex = plotx < 360 & 0 < plotx;
                    plotx = plotx(dex);
                    ploty = ploty(dex);
                    scattercolors = scattercolors(dex,:);
                    %set(h2(j), {'color'}, linecolors);
                    %scatter(ax1,reshape(plotx, [numel(plotx) 1])-poloidalextent*k, reshape(ploty, [numel(ploty) 1]), 5, reshape(plotcolors, [size(plotcolors,1) 3]));
                    
                    plot(ax1,plotx,ploty, 'Color', linecolors(j,:));
                    scatter(ax1,plotx,ploty, 5, scattercolors);
                    
                    if pspace
                        plot(ax2,pll(i-tail:i,j)./amu,pperp(i-tail:i,j)./amu, 'LineWidth', 5, 'Color', linecolors(j,:));
                        for m = plot_torrevolutions(j,1):plot_torrevolutions(j,1) %Frames for toroidal crossings
                            ploty = phi(i-tail:i,j)-toroidalextent*m;
                            ploty = ploty(dex);
                            %ploty = mod(ploty(dex),toroidalextent); %Mod should be removed by previous calculations
                            plot(ax3, plotx, ploty, 'Color', linecolors(j,:));
                        end
                    end
                    
                    %scatter(ax1,plotx(:,j)-poloidalextent*plotrevolutions(j), ploty(:,j), 5, reshape(plotcolors, [size(plotcolors,1) 3]));
                end
            end
            
            %end
            %scatter(ax1,reshape(plotx, [numel(plotx) 1]), reshape(ploty, [numel(ploty) 1]), 5, reshape(plotcolors, [size(plotcolors,1) 3]));
            %                 radii = .001 + mod(phi(i,:),360)/20000;
            %                 d = 2 .* radii;
            %                 da = daspect(ax1);
            %                 px = mean(plotx(end-10:end,:),1) - radii.*da(1)/2;
            %                 py = mean(ploty(end-10:end,:),1) - radii.*da(2);
            %
            %                 for k = 1:size(plotx,2)
            %                 r(k) = rectangle(ax1,'Position',[px(k) py(k) d(k)*da(1)/2 d(k)*da(2)],'Curvature',[1,1]);
            %                 end
            if pspace
                %h3 = plot(ax2,pll(i-tail:i,:)./amu,pperp(i-tail:i,:)./amu, 'LineWidth', 5);
                %set(h3, {'color'}, linecolors);
                ax2.XLim = [-3E6 3E6];
                ax2.YLim = [0 4E6];
                xlabel(ax2,'Parallel Mom. [u*m/s]');
                ylabel(ax2,'Perp. Mom. [u*m/s]')
                %title(ax2,'Phase space');
                set(ax2,'FontSize',20);
                %h5 = plot(ax3, plotx(:,:), mod(phi(i-tail:i,:),toroidalextent));
                ax3.XLim = [0 poloidalextent];
                ax3.YLim = [0 toroidalextent];
                xlabel(ax3, 'Poloidal angle \theta [deg]');
                ylabel(ax3, 'Toroidal angle \phi [deg]');
                set(ax3,'FontSize',20);
                %set(h5, {'color'}, linecolors);
            end
            ax1.YLim = [0 1.5];
            ax1.XLim = [0 poloidalextent];
            ax1.CLim = [min(b, [], 'all'), max(b, [], 'all')];
            c = colorbar(ax1,'south');
            c.Label.String = 'B Field [T]';
            %scatter(ax1, deg2rad(theta(i,:),rho(i,:));
            set(ax1,'FontSize',24);
            
            xlabel(ax1, 'Poloidal angle \theta [deg]');
            ylabel(ax1, 'Flux label \rho [a.u.]');
            %title(ax1,'Particle Orbits');
            frame = getframe(fig); %for writing
            writeVideo(v,frame);
            if lmovie
            cla(ax1)
            if pspace
                cla(ax2)
                cla(ax3)
            end
            end
            
        end
        %         else
        %             plot(ax1,theta,rho);
        %             set(gca,'FontSize',24);
        %             title('Particle Orbits');
        %             ax1.YLim = [0 1.5];
        %             ax1.XLim = [0 poloidalextent];
        %             ax1.CLim = [min(b, [], 'all'), max(b, [], 'all')];
        %             c = colorbar(ax1,'south');
        %             c.Label.String = 'B Field [T]';
        %             %scatter(ax1, deg2rad(theta(i,:),rho(i,:));
        %             set(ax1,'FontSize',24);
        %
        %             xlabel(ax1, 'Poloidal angle \theta [deg]');
        %             ylabel(ax1, 'Flux label \rho [a.u.]');
        %             title(ax1,'Particle Orbits');
        %         end
    case {'vll','pll'}
        
        % Make plot
        
        plot(ax1,pll./amu,pperp./amu);
        set(gca,'FontSize',24);
        xlabel('Parallel Momentum [u*m/s]');
        ylabel('Perpendicular Momentum [u*m/s]')
        title('Phase space');
    case 'rz'
        % Make plot
        plot(ax1, r,z);
        set(gca,'FontSize',24);
        xlabel('R [m]');
        ylabel('Z [m]')
        title('R-Z Overview');
        
end
end

