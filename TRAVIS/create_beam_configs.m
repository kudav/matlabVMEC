function beamConfigs = create_beam_configs(startRPZ, endRPZ, numSteeringAngles, angleStep)
% CREATE_BEAM_CONFIGS Creates beam configurations for steering angles
%
% Inputs:
%   startRPZ - matrix, start positions [r, z, phi] 
%   endRPZ - matrix, end positions [r, z, phi]
%   numSteeringAngles - number of steering angles to generate
%   angleStep - angular step between beams in degrees
%
% Returns:
%   beamConfigs - structure array with beam configurations

    if nargin < 4
        numSteeringAngles = 1;
        angleStep = 0;
    end
    
    numBeams = size(startRPZ, 1);
    beamConfigs = repmat(struct('start', [], 'end', []), numBeams * numSteeringAngles, 1);
    
    idx = 1;
    for i = 1:numBeams
        for j = 0:numSteeringAngles-1
            beamConfigs(idx).start = startRPZ(i, :);
            beamConfigs(idx).end = endRPZ(i, :) + [0, 0, j * angleStep];
            idx = idx + 1;
        end
    end
end