function generate_steered_beams(inputFilename, outputFilename, startRPZ, endRPZ, dist, refaxis, angle, n, m)
% GENERATE_STEERED_BEAMS Generates multiple beams with different steering angles
%
% Inputs:
%   inputFilename - string, path to input XML file
%   outputFilename - string, path to output XML file
%   startRPZ - matrix, start positions [r, z, phi] for each beam
%   endRPZ - matrix, end positions [r, z, phi] for each beam  
%   dist - scalar, distance parameter
%   refaxis - vector, reference axis for steering (default [0, 0, 1])
%   angle - scalar, steering angle in degrees
%   n - scalar, number of beams in radial direction
%   m - scalar, number of beams in toroidal direction
%
% Example:
%   startRPZ = [5.44, 34.69, 0.78; 5.50, 34.70, 0.80];
%   endRPZ = [5.06, 31.70, -0.10; 5.10, 31.75, -0.05];
%   generate_steered_beams('input.xml', 'output.xml', startRPZ, endRPZ, 0.12, [0, 0, 1], 1.5, 1, 1);

    % Validate inputs
    if nargin < 9
        error('Missing required arguments');
    end
    
    if size(startRPZ, 1) ~= size(endRPZ, 1)
        error('startRPZ and endRPZ must have the same number of rows');
    end
    
    % Read the input XML structure
    input = readstruct(inputFilename);
    
    % Get number of original beams and create new ones
    numbeams = size(startRPZ, 1);
    
    % Store original beam for copying
    originalBeam = input.ECRHsystem.Beam;
    
    % Generate additional beams
    for i = 1:numbeams
        % Copy the original beam
        tempBeam = originalBeam;
        
        % Add new beam to the structure
        if isfield(input.ECRHsystem, 'Beam')
            if isnumeric(input.ECRHsystem.Beam)
                % If Beam is numeric, convert to cell array
                input.ECRHsystem.Beam = {input.ECRHsystem.Beam};
            end
            
            % Add new beam at the end
            input.ECRHsystem.Beam{end+1} = tempBeam;
            
            % Update beam properties
            beamIndex = end;
            input.ECRHsystem.Beam{beamIndex}.origin.Text = sprintf('%.2f %.2f %.2f', startRPZ(i,:));
            input.ECRHsystem.Beam{beamIndex}.direction.Text = sprintf('%.2f %.2f %.2f', endRPZ(i,:));
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