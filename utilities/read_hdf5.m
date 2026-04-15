function data = read_hdf5(filename, subgroups)
%READ_HDF5 Returns the contents of an HDF5 file as a structure
%   data = read_hdf5(filename)
%   data = read_hdf5(filename, subgroups)
%
%   Optional input:
%     subgroups : names/paths of subgroups to read. If supplied, ONLY those
%                 groups are read (root datasets are NOT read).
%                 Examples:
%                   read_hdf5('file.h5', {'/group1','/group2/sub'})
%                   read_hdf5('file.h5', {'group1','group2/sub'})  % leading '/' optional

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

data = struct();

% Parse optional subgroup list
readOnlySomeGroups = (nargin >= 2) && ~isempty(subgroups);
if readOnlySomeGroups
    if ischar(subgroups) || (isstring(subgroups) && isscalar(subgroups))
        subgroups = cellstr(subgroups);
    elseif isstring(subgroups)
        subgroups = cellstr(subgroups);
    end
    % normalize: ensure leading '/' and remove trailing '/'
    for k = 1:numel(subgroups)
        g = subgroups{k};
        if isempty(g), continue; end
        if g(1) ~= '/', g = ['/' g]; end
        if numel(g) > 1 && g(end) == '/', g(end) = []; end
        subgroups{k} = g;
    end
end

ngroups = length(data_info.Groups);
nvars   = length(data_info.Datasets);

% Read root datasets ONLY if no subgroup filter was provided
if ~readOnlySomeGroups
    for i = 1:nvars
        name_local = data_info.Datasets(i).Name;
        name_local = strrep(name_local,' ','_');
        data.(name_local) = h5read(filename, ['/' data_info.Datasets(i).Name]);

        natts = length(data_info.Datasets(i).Attributes);
        for j = 1:natts
            att_name_local = data_info.Datasets(i).Attributes(j).Name;
            if ~contains(att_name_local, name_local)
                att_name_local = [name_local '_att_' strrep(att_name_local,' ','_')];
            end
            data.(att_name_local) = data_info.Datasets(i).Attributes(j).Value;
        end
    end
end

% Read groups
if readOnlySomeGroups
    % Read only requested groups (by path)
    for k = 1:numel(subgroups)
        grpPath = subgroups{k};
        if isempty(grpPath) || strcmp(grpPath,'/'), continue; end
        try
            field_name = strrep(grpPath(2:end), ':', '_'); % drop leading '/'
            field_name = matlab.lang.makeValidName(field_name);
            data.(field_name) = getGroup(filename, grpPath);
        catch ME
            warning('read_hdf5:GroupReadFailed', ...
                'Failed to read group "%s": %s', grpPath, ME.message);
        end
    end
else
    % Original behavior: read all root groups
    if ngroups > 0
        for i = 1:ngroups
            group_name = data_info.Groups(i).Name;
            group_name = strrep(group_name,':','_');
            data.(group_name(2:end)) = getGroup(filename, data_info.Groups(i).Name);
        end
    end
end
end

function data = getGroup(filename, rootdir)
data_info = h5info(filename, rootdir);
data = walkGroup(filename, data_info);
% ngroups   = length(data_info.Groups);
% nvars     = length(data_info.Datasets);
% 
% data = struct();
% 
% % Datasets in this group
% for i = 1:nvars
%     name_local = data_info.Datasets(i).Name;
%     name_local = strrep(name_local,' ','_');
% 
%     data.(data_info.Datasets(i).Name) = h5read(filename, [rootdir '/' data_info.Datasets(i).Name]);
% 
%     natts = length(data_info.Datasets(i).Attributes);
%     for j = 1:natts
%         att_name_local = data_info.Datasets(i).Attributes(j).Name;
%         if ~startsWith(att_name_local, name_local)
%             att_name_local = [name_local '_att_' strrep(att_name_local,' ','_')];
%         end
%         data.(att_name_local) = data_info.Datasets(i).Attributes(j).Value;
%     end
% end
% 
% % Subgroups (recursive)
% for i = 1:ngroups
%     group_name = data_info.Groups(i).Name;
%     dex = strfind(group_name,'/');
%     group_name = group_name(dex(end)+1:end);
%     try
%         data.(group_name) = getGroup(filename, ['/' data_info.Groups(i).Name]);
%     catch
%         group_name = ['GID_' group_name];
%         data.(group_name) = getGroup(filename, ['/' data_info.Groups(i).Name]);
%     end
% end
end

function data = walkGroup(filename, info)
data = struct();
for i = 1:length(info.Datasets)
    name_local = strrep(info.Datasets(i).Name,' ','_');
    fullpath   = [info.Name '/' info.Datasets(i).Name];
    fullpath   = strrep(fullpath,'//','/');
    data.(name_local) = h5read(filename, fullpath);
    for j = 1:length(info.Datasets(i).Attributes)
        att = info.Datasets(i).Attributes(j).Name;
        if ~contains(att, name_local)
            att = [name_local '_att_' strrep(att,' ','_')];
        end
        data.(att) = info.Datasets(i).Attributes(j).Value;
    end
catch
    % Fall back to path-based dedup if H5O is unavailable
    key = path;
end
for i = 1:length(info.Groups)
    gname = info.Groups(i).Name;
    gname = gname(find(gname=='/',1,'last')+1:end);
    gname = strrep(gname,':','_');
    try
        data.(gname) = walkGroup(filename, info.Groups(i));
    catch
        data.(['GID_' gname]) = walkGroup(filename, info.Groups(i));
    end
end
end