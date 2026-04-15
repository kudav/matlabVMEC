function data = read_hdf5(filename)
%READ_HDF5 Returns the contents of an HDF5 file as a structure
%   The READ_HDF5 function reads an HDF5 file and returns the contents of
%   that file as the fields of a structure.  Groups are treated as elements
%   of their parent structure.  If the file file cannot be opened a -1 is
%   returned.
%
%   Example
%       data=read_hdf5('input.h5');
%
%   Version 1.3
%   Maintained by: Samuel Lazerson (lazerson@pppl.gov)
%   Date  05/02/2012


% Try to read the file first
try
    data_info = h5info(filename,'/');
catch h5info_error
    data=-1;
    disp(['ERROR: Opening HDF5 File: ' filename]);
    disp(['  -identifier: ' h5info_error.identifier]);
    disp(['  -message:    ' h5info_error.message]);
    disp('      For information type:  help read_hdf5');
    return
end

% h5info('/') already returns the entire tree recursively.
% Walk that cached tree once instead of re-calling h5info per subgroup,
% which on deeply nested files caused effectively-exponential re-scans.
data = walkGroup(filename, data_info);
return
end

function data = walkGroup(filename, info)
data = struct();

for i = 1:length(info.Datasets)
    name_local = strrep(info.Datasets(i).Name,' ','_');
    fullpath   = [info.Name '/' info.Datasets(i).Name];
    fullpath   = strrep(fullpath,'//','/');
    data.(name_local) = h5read(filename, fullpath);
    natts = length(info.Datasets(i).Attributes);
    for j = 1:natts
        att_name_local = info.Datasets(i).Attributes(j).Name;
        if ~contains(att_name_local,name_local)
            att_name_local = [name_local '_att_' strrep(att_name_local,' ','_')];
        end
        data.(att_name_local) = info.Datasets(i).Attributes(j).Value;
    end
end

for i = 1:length(info.Groups)
    group_name = info.Groups(i).Name;
    dex = strfind(group_name,'/');
    group_name = group_name(dex(end)+1:end);
    group_name = strrep(group_name,':','_');
    try
        data.(group_name) = walkGroup(filename, info.Groups(i));
    catch
        group_name = ['GID_' group_name];
        data.(group_name) = walkGroup(filename, info.Groups(i));
    end
end
return
end
