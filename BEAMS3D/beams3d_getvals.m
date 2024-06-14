function values = beams3d_getvals(data, r, phi, z, varargin)
    %BEAMS3D_GETVALS Interpolates 3D fields to the queried r, phi, z coordinates
    %   This function takes in r, phi, z coordinates and a list of field names,
    %   and interpolates the 3D fields corresponding to these names to the 
    %   queried (r, phi, z) coordinates.
    %   Valid field names are: 'TE', 'NE', 'TI', 'NI', 'S_ARR', 'U_ARR'.

    % Initialize output structure
    values = struct();

    % List of valid field names
    valid_field_names = {'TE', 'NE', 'TI', 'NI', 'S_ARR', 'U_ARR'};
    
    % Iterate over each requested field name in varargin
    for i = 1:length(varargin)
        field_name = varargin{i};
        
        % Check if the field name is valid
        if ismember(field_name, valid_field_names)
            % Interpolate the field
            %field_data = getfield(data, field_name);  % Access field data
            interpolated_values = interp3(data.raxis, data.phiaxis, data.zaxis, ...
                                          permute(data.(field_name), [2 1 3]), ...
                                          r, mod(phi, data.phiaxis(end)), z, 'linear', 0);
            % Store the interpolated values in the output structure
            values.(field_name) = interpolated_values;
        else
            error('Invalid field name: %s', field_name);
        end
    end
end
