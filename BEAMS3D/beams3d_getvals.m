function values = beams3d_getvals(data, r, phi, z, varargin)
%BEAMS3D_GETVALS Interpolates 3D fields to the queried r, phi, z coordinates
%   This function takes in r, phi, z coordinates and a list of field names,
%   and interpolates the 3D fields corresponding to these names to the
%   queried (r, phi, z) coordinates.
%   Valid field names are: 'TE', 'NE', 'TI', 'NI', 'S_ARR', 'U_ARR'.

if ischar(data)
    b3d_filename=data;
    data={};
    data.raxis = h5read(b3d_filename,'/raxis');
    data.phiaxis = h5read(b3d_filename,'/phiaxis');
    data.zaxis = h5read(b3d_filename,'/zaxis');
end
% Initialize output structure
values = struct();

% List of valid field names
valid_field_names = {'TE', 'NE', 'TI',...
    'S_ARR', 'U_ARR','POT_ARR', 'ZEFF_ARR',...
    'B_R', 'B_PHI','B_Z','MODB'};

% Iterate over each requested field name in varargin
for i = 1:length(varargin)
    field_name = varargin{i};
    % Check if the field name is valid
    if ismember(field_name, valid_field_names)
        if ~isfield(data,field_name)%Load field if not available
            if strcmp(field_name,'MODB')
                data.(field_name)=sqrt(data.B_R.^2+data.B_PHI.^2+data.B_Z.^2);
            else
                data.(field_name)= h5read(b3d_filename,['/',field_name]);
            end
        end
        % Interpolate the field
        %field_data = getfield(data, field_name);  % Access field data
        interpolated_values = interp3(data.raxis, data.phiaxis, data.zaxis, ...
            permute(data.(field_name), [2 1 3]), ...
            r, mod(phi, data.phiaxis(end)), z, 'cubic', 0);
        % Store the interpolated values in the output structure
        values.(field_name) = interpolated_values;
    else
        error('Invalid field name: %s', field_name);
    end
end
end
