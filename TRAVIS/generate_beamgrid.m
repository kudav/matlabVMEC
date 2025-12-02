function generate_beamgrid(filename, n, m, dist, angle,refaxis)
if isempty(filename)
    filename='input.xml';
end
if isempty(refaxis)
    refaxis = [0, 0, 1]; % Default reference axis if not provided
end
input = readstruct(filename);
%%
numbeams=numel(input.ECRHsystem.Beam);
k=numbeams-1;
for i=1:numbeams
    tempBeam=input.ECRHsystem.Beam(i);
    startRPZ = sscanf(tempBeam.origin.Text, '%f', 3);%[.61 deg2rad(31.34) 0.83];
    endRPZ   =sscanf(tempBeam.direction.Text, '%f', 3);%[5.16 deg2rad(28.70) -0.02];
    startRPZ(2)=deg2rad(startRPZ(2));
    endRPZ(2)=deg2rad(endRPZ(2));
    [startGrid, endGrid] = parallelBeamGrid(startRPZ, endRPZ, 2, 3, dist, 'RefAxis',refaxis,'FocusAngleDeg',angle);%,'FocusAngleDeg', 5
    startGrid(:,:,2)=rad2deg(startGrid(:,:,2));
    endGrid(:,:,2)=rad2deg(endGrid(:,:,2));
    for j=1:n*m
        k=k+1;
        [h,g]=ind2sub([n,m],j);
        input.ECRHsystem.Beam(end+1)=tempBeam;
        input.ECRHsystem.Beam(end).origin.Text = sprintf("%.2f %.2f %.2f", startGrid(g,h,:));
        input.ECRHsystem.Beam(end).direction.Text = sprintf("%.2f %.2f %.2f", endGrid(g,h,:));
        input.ECRHsystem.Beam(end).idAttribute=k;
        input.ECRHsystem.Beam(end).enabledAttribute=1;
        input.ECRHsystem.Beam(end).nameAttribute=strcat(input.ECRHsystem.Beam(end).nameAttribute, " - Grid 1");
        input.ECRHsystem.Beam(end).power.Text=tempBeam.power.Text/(n*m);
        input.ECRHsystem.Beam(end).nCircles.Text=0;
    end
    input.ECRHsystem.Beam(i).enabledAttribute=0;
end

outputFilename = strrep(filename, '.xml', '_grid.xml');
% Write the modified structure to the output XML file
writestruct(input, outputFilename,'StructNodeName','ECRHheating');

% Remove all <Text> and </Text> strings from the output XML file while keeping the text between them unaffected
outputContent = fileread(outputFilename);
outputContent = regexprep(outputContent, '\s*<Text>\s*(.*?)\s*<\/Text>\s*', '$1');
fid = fopen(outputFilename, 'w');
fwrite(fid, outputContent);
fclose(fid);
end