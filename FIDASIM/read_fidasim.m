function fidasim_out = read_fidasim(runid, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if (strcmp(runid(end-1:end),'h5'))
    disp('ERROR: Only give runid (dist and eq runids are loaded automatically!');
    disp(['       runidname: ' runid]);
    return
end

ldist=1;
lgeom=1;
lspec=1;
lneut=1;
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