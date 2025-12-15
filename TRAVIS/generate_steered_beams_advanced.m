function generate_steered_beams_advanced(inputFilename, outputFilename, beamConfigs, dist, refaxis, angle, n, m)
% GENERATE_STEERED_BEAMS_ADVANCED Advanced version with beam configuration structure
%
% Inputs:
%   inputFilename - string, path to input XML file
%   outputFilename - string, path to output XML file
%   beamConfigs - structure array or cell array with beam configurations
%   dist - scalar, distance parameter
%   refaxis - vector, reference axis for steering (default [0, 0, 1])
%   angle - scalar, steering angle in degrees
%   n - scalar, number of beams in radial direction
%   m - scalar, number of beams in toroidal direction
%
% Example:
%   beamConfigs = struct('start', [5.44, 34.69, 0.78], 'end', [5.06, 31.70, -0.10]);
%   generate_steered_beams_advanced('input.xml', 'output.xml', beamConfigs, 0.12, [0, 0, 1], 1.5, 1, 1);

    % Validate inputs
    if nargin < 8
        error('Missing required arguments');
    end
    
    % Read the input XML structure
    input = readstruct(inputFilename);
    
    % Get number of beam configurations
    if isstruct(beamConfigs)
        numbeams = length(beamConfigs);
    elseif iscell(beamConfigs)
        numbeams = length(beamConfigs);
    else
        numbeams = size(beamConfigs, 1);
    end
    
    % Store original beam for copying
    originalBeam = input.ECRHsystem.Beam;
    
    % Generate additional beams based on configurations
    for i = 1:numbeams
        % Copy the original beam
        tempBeam = originalBeam;
        
        % Determine beam configuration
        if isstruct(beamConfigs)
            startRPZ = beamConfigs(i).start;
            endRPZ = beamConfigs(i).end;
        elseif iscell(beamConfigs)
            startRPZ = beamConfigs{i}.start;
            endRPZ = beamConfigs{i}.end;
        else
            startRPZ = beamConfigs(i, :);
            endRPZ = beamConfigs(i + numbeams, :); % Assuming second half contains end positions
        end
        
        % Add new beam to the structure
        if isfield(input.ECRHsystem, 'Beam')
            if isnumeric(input.ECRHsystem.Beam)
                % If Beam is numeric, convert to cell array
                input.ECRHsystem.Beam = {input.ECRHsystem.Beam};
            end
            
            % Add new beam at the end
            input.ECRHsystem.Beam{end+1} = tempBeam;
            
            % Update beam properties
            %beamIndex = end;
            input.ECRHsystem.Beam{beamIndex}.origin.Text = sprintf('%.2f %.2f %.2f', startRPZ);
            input.ECRHsystem.Beam{beamIndex}.direction.Text = sprintf('%.2f %.2f %.2f', endRPZ);
            input.ECRHsystem.Beam{beamIndex}.idAttribute = i + 2;
            input.ECRHsystem.Beam{beamIndex}.enabledAttribute = 1;
            input.ECRHsystem.Beam{beamIndex}.nameAttribute = sprintf('BEAM %d', i);
            
            % Adjust power for multiple beams
            if isfield(tempBeam, 'power') && isfield(tempBeam.power, 'Text')
                input.ECRHsystem.Beam{beamIndex}.power.Text = num2str(str2double(tempBeam.power.Text) / (n * m));
            end
            
            % Set other beam parameters
            if isfield(tempBeam, 'nCircles') && isfield(tempBeam.nCircles, 'Text')
                input.ECRHsystem.Beam{beamIndex}.nCircles.Text = '0';
            end
        end
    end
    
    % Disable original beams (first two beams)
    if isfield(input.ECRHsystem, 'Beam') && isnumeric(input.ECRHsystem.Beam)
        input.ECRHsystem.Beam(1).enabledAttribute = 0;
        input.ECRHsystem.Beam(2).enabledAttribute = 0;
    elseif isfield(input.ECRHsystem, 'Beam') && iscell(input.ECRHsystem.Beam)
        if length(input.ECRHsystem.Beam) >= 2
            input.ECRHsystem.Beam{1}.enabledAttribute = 0;
            input.ECRHsystem.Beam{2}.enabledAttribute = 0;
        end
    end
    
    % Write the modified structure to the output XML file
    writestruct(input, outputFilename, 'StructNodeName', 'ECRHheating');
    
    % Remove all <Text> and </Text> tags while preserving content
    outputContent = fileread(outputFilename);
    outputContent = regexprep(outputContent, '<Text>(.*?)</Text>', '$1');
    outputContent = regexprep(outputContent, '\s*<Text>\s*(.*?)\s*</Text>\s*', '$1');
    
    % Write back to file without Text tags
    fid = fopen(outputFilename, 'w');
    if fid == -1
        error('Cannot open output file: %s', outputFilename);
    end
    fwrite(fid, outputContent, 'char');
    fclose(fid);
    
    fprintf('Successfully generated steered beams. Output written to %s\n', outputFilename);
end
