function fidasim_out = read_fidasim(runid, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if (strcmp(runid(end-1:end),'h5'))
    disp('ERROR: Only give runid (dist and eq runids are loaded automatically!');
    disp(['       runidname: ' runid]);
    return
end
b3d_S=[];
if nargin > 1
    i = 1;
    ldist=0;
    lgeom=0;
    lspec=0;
    lneut=0;
    leq=0;
    lweight=0;
    lbirth=0;
    while i < nargin
        switch lower(varargin{i})
            case 'birth'
                lbirth=1;
            case 'dist'
                ldist=1;
            case 'geom'
                lgeom=1;
            case 'spec'
                lspec=1;
            case 'neut'
                lneut=1;
            case 'eq'
                leq=1;
            case 'weight'
                lweight=1;
            case 'b3d'
                if nargin==3
                    ldist=1;
                    lgeom=1;
                    lspec=1;
                    lneut=1;
                    leq=1;
                    lweight=1;
                    lbirth=1;    
                end
                i = i+1;
                b3d_filename=varargin{i};
                b3d_S= h5read(b3d_filename,'/S_ARR');
                r = h5read(b3d_filename,'/raxis');
                phi = h5read(b3d_filename,'/phiaxis');
                z = h5read(b3d_filename,'/zaxis');
        end
        i=i+1;
    end
else
    ldist=1;
    lgeom=1;
    lspec=1;
    lneut=1;
    leq=1;
    lweight=1;
    lbirth=1;
end

nml_name=[runid,'_inputs.dat'];
if isfile(nml_name)
    input=read_namelist(nml_name,'fidasim_inputs');
    % Get all field names in the structure
    fields = fieldnames(input);

    % Loop through each field
    for i = 1:numel(fields)
        % Get the current field's value
        currentValue = input.(fields{i});

        % Check if the field contains a character array
        if ischar(currentValue)
            % Replace the old substring with the new one
            input.(fields{i}) = strrep(currentValue, '%%FIDASIM_PATH%%NUMBER', runid(1:2));
        end
    end
    runid=input.runid;
    dist_name=splitfn(input.distribution_file);
    eq_name=splitfn(input.equilibrium_file);
    geom_name=splitfn(input.geometry_file);
    neut_name=splitfn(input.neutrals_file);
    fidasim_out.input=input;
else
    disp('Could not find input namelist, using defaults!')
    dist_name = [runid, '_distribution.h5'];
    eq_name = [runid,'_equilibrium.h5'];
    neut_name = [runid, '_neutrals.h5'];
    geom_name = [runid,'_geometry.h5'];
end

birth_name = [runid,'_birth.h5'];
weight_name = [runid,'_fida_weights.h5'];
spec_name = [runid,'_spectra.h5'];

if isfile(dist_name)&&ldist
    disp([' Reading file: ' dist_name]);
    dist = read_hdf5(dist_name);
    fidasim_out.dist=dist;
elseif ~ldist
    disp('Skipping distribution')
else
    disp('ERROR: Distribution file not found, check runidname!');
end

if isfile(eq_name)&&leq
    disp([' Reading file: ' eq_name]);
    eq = read_hdf5(eq_name);
    if ~isempty(b3d_S)
        [rgrid,phigrid,zgrid]=ndgrid(eq.fields.r/100,eq.fields.phi,eq.fields.z/100);
        eq.fields.s = interp3(r,phi,z,...
            permute(sqrt(b3d_S),[2 1 3]),rgrid,mod(phigrid,phi(end)),zgrid,'linear',NaN);
        eq.fields.s=permute(eq.fields.s,[1 3 2]);
    end    
    fidasim_out.eq=eq;

elseif ~leq
    disp('Skipping equilibrium')
else
    disp('ERROR: Equilbirium file not found!');
end

if isfile(weight_name)&&lweight
    disp([' Reading file: ' weight_name]);
    weight= read_hdf5(weight_name);
    fidasim_out.weight=weight;
elseif ~lweight
    disp('Skipping weight fcn.')
else
    disp('ERROR: Weight file not found!');
end

if isfile(neut_name)&&lneut
    disp([' Reading file: ' neut_name]);
    neut = read_hdf5(neut_name);
    fidasim_out.neut=neut;
elseif ~lneut
    disp('Skipping neutrals')
else
    disp('ERROR: Neutrals file not found!');
end
if isfile(spec_name)&&lspec
    disp([' Reading file: ' spec_name]);
    spec = read_hdf5(spec_name);
    fidasim_out.spec=spec;
elseif ~lspec
    disp('Skipping spectra')
else
    disp('ERROR: Spectra file not found!');
end

if isfile(geom_name)&&lgeom
    disp([' Reading file: ' geom_name]);
    geom = read_hdf5(geom_name);
    if isfield(fidasim_out.eq.fields,'s')
        disp('Calculating rho values of LOS')
         % Number of lines of sight in spec
        nchan = geom.spec.nchan;
        
        % Initialize array to hold closest points
        closest_points = zeros(3, nchan);
        
        % Loop through each line of sight in spec
        for i = 1:nchan
            % Get the line of sight direction and a point on the line
            % spec_line_dir = fidasim_out.geom.spec.axis(:,i);
            % spec_point = fidasim_out.geom.spec.lens(:,i);         
            % Find the closest point on nbi_axis to the current line of sight
            closest_points(:,i) = closest_point_on_line(geom.nbi.axis,geom.nbi.src, geom.spec.axis(:,i), geom.spec.lens(:,i));
        end
        geom.spec.closest_points=closest_points;
        [closest_points(2,:),closest_points(1,:),closest_points(3,:)]=cart2pol(closest_points(1,:),closest_points(2,:),closest_points(3,:));
        geom.spec.closest_points_cyl=closest_points;
        geom.spec.rho = interp3(r*100,phi,z*100,...
            permute(sqrt(b3d_S),[2 1 3]),closest_points(1,:),mod(closest_points(2,:),phi(end)),closest_points(3,:),'linear',NaN);
        %geom.spec.rho = interp3(eq.fields.r,eq.fields.z,eq.fields.phi,...
        %    permute(eq.fields.s,[2 1 3]),closest_points(1,:),closest_points(3,:),mod(closest_points(2,:),phi(end)),'linear');       
    end
    fidasim_out.geom=geom;
elseif ~lgeom
    disp('Skipping geometry')
else
    disp('ERROR: Geometry file not found!');
end
if isfile(birth_name)&&lbirth
    disp([' Reading file: ' birth_name]);
    birth = read_hdf5(birth_name);
    fidasim_out.birth=birth;
elseif ~lbirth
    disp('Skipping births')
else
    disp('ERROR: Birth file not found!');
end
end

function end_el=splitfn(stringvar)
splitResult = split(stringvar, '/');
end_el = splitResult{end};
end



function closest_point = closest_point_on_line(axis_dir, axis_point, line_dir, line_point)
    % axis_dir: 3x1 vector representing the direction of the nbi axis
    % axis_point: 3x1 vector representing a point on the nbi axis
    % line_dir: 3x1 vector representing the direction of the line of sight
    % line_point: 3x1 vector representing a point on the line of sight
    
    % Calculate the vector between the two points
    w0 = axis_point - line_point;
    
    % Calculate the coefficients for the lines
    a = dot(axis_dir, axis_dir);
    b = dot(axis_dir, line_dir);
    c = dot(line_dir, line_dir);
    d = dot(axis_dir, w0);
    e = dot(line_dir, w0);
    
    % Calculate the parameters that minimize the distance
    sc = (b*e - c*d) / (a*c - b^2);
    tc = (a*e - b*d) / (a*c - b^2);
    
    % Calculate the closest points on both lines
    closest_point_axis = axis_point + sc * axis_dir;
    closest_point_line = line_point + tc * line_dir;
    
    % The closest point on the line of sight to the axis line
    closest_point = closest_point_line;
end