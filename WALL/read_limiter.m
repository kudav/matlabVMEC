function data = read_limiter( filename )
%READ_LIMITER Read's a limiting surface file
%   READ_LIMITER reads a limiter data file.  The limiter data file is a
%   file with the format:
%
%MACHINE:  W-7X_limiter
%DATE:  07-01-14
%07612 14399
%    5.5057190000E+00    -7.4592000000E-02    -4.4169700000E-01
%    5.5057190000E+00     7.4592000000E-02     4.4169700000E-01
%    5.5063620000E+00    -7.2601000000E-02    -4.4126700000E-01
%    5.5063620000E+00     7.2601000000E-02     4.4126700000E-01
%    5.5165410000E+00    -7.5732000000E-02    -4.1455700000E-01
%    5.5165410000E+00     7.5732000000E-02     4.1455700000E-01
%         .                      .                     .
%         .                      .                     .
%         .                      .                     .
%    1 2 3
%    4 5 6
%    . . .
%
%   Here the string after MACHINE is the name of the machine.  The DATE
%   string contains information about the date of the file.  The next line
%   indicates the number of verticies and number of triangular cells.  What
%   follows are a list of vertices (x,y,z) then the list of points
%   composing each cell.
%
%   Example:
%       lim_data=read_limiter('limiter_trimesh.dat');
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           1/10/13
% Check arguments
if nargin<1
    disp('ERROR: read_vessel requires filename');
    return
end
% Open File
fid=fopen(filename,'r');
% Read Header
header_line1=fgetl(fid);
header_line2=fgetl(fid);
data.machine=strtrim(header_line1(strfind(header_line1,':')+1:numel(header_line1)));
data.date=strtrim(header_line2(strfind(header_line2,':')+1:numel(header_line2)));
temp=fscanf(fid,'%d',2);
data.nvertex = temp(1);
data.nfaces = temp(2);
% Read dataset
data.coords=fscanf(fid,'%E %E %E',[3 data.nvertex]);
data.faces=fscanf(fid,'%d %d %d',[3 data.nfaces]);
fclose(fid);
data.datatype='limiter_trimesh';


end

