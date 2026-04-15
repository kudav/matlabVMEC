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
%   Version 1.4
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

% h5info('/') already returns the entire tree recursively, so we walk that
% cached struct rather than re-calling h5info per subgroup. We also keep a
% visited-set keyed on the canonical object address (via H5O.get_info), so
% that hard links / soft links / cycles do not cause the same group to be
% read repeatedly.
fid     = H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
visited = containers.Map('KeyType','char','ValueType','logical');
cleanupObj = onCleanup(@() H5F.close(fid));

data = walkGroup(filename, fid, data_info, visited);
return
end

function data = walkGroup(filename, fid, info, visited)
data = struct();

key = objKey(fid, info.Name);
if ~isempty(key)
    visited(key) = true;
end

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
    gpath  = info.Groups(i).Name;
    subkey = objKey(fid, gpath);
    if ~isempty(subkey) && isKey(visited, subkey)
        continue   % alias / cycle — already read
    end
    group_name = gpath(find(gpath=='/',1,'last')+1:end);
    group_name = strrep(group_name,':','_');
    try
        data.(group_name) = walkGroup(filename, fid, info.Groups(i), visited);
    catch
        group_name = ['GID_' group_name];
        data.(group_name) = walkGroup(filename, fid, info.Groups(i), visited);
    end
end
return
end

function key = objKey(fid, path)
% Return a string identifying the underlying HDF5 object so that aliases
% (multiple hard links to one group) collapse to the same key.
key = '';
try
    oinfo = H5O.get_info_by_name(fid, path, 'H5P_DEFAULT');
    if isstruct(oinfo) && isfield(oinfo,'addr')
        key = sprintf('%lu', uint64(oinfo.addr));
    elseif isstruct(oinfo) && isfield(oinfo,'token')
        key = sprintf('%d_', oinfo.token(:));
    end
catch
    % Fall back to path-based dedup if H5O is unavailable
    key = path;
end
return
end
