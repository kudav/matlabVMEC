function data = read_travis
%read_travis Reads the files in a TRAVIS output directory
%   The READ_TRAVIS routine reads the files in the current directory into a
%   a structure for plotting.  For the beamtrace cell arrays are created
%   which are [nrays,nbeams] in size (nbeams being the number of beamtrace
%   files in the directory).
%
%   Example:
%       travis_data=read_travis();
%
%   Written by:     S.Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:        1.0
%   Date:           2/5/23


%files=get_files_sorted('beamtrace*');
[~,files] = get_files_sorted('beamtrace*');

columnames = {'Nray','path','x','y','z','nx','ny','nz','rho','ne','te','B', ...
    'Nper','Npar','Nperc','Nparc','damp0','damp','tau0','tau', ...
    'wray0','wray','effcd','dPdlp','dPdlt','dIdlp','dIdlt','cosVk', ...
    'wReEx','wImEx','wReEy','wImEy','wReEz','wImEz', ...
    'cReEx','cImEx','cReEy','cImEy','cReEz','cImEz', ...
    'uminP','ucenP','umaxP','uminj','ucenj','umaxj','upar1','upar2', ...
    'Bx','By','Bz','Npass'};

for i = 1:numel(files)
    S = importdata(files(i).name, ' ', 1);
    temp = S.data;
    nrays = max(temp(:,1));
    %nsteps=accumarray(temp(:,1),1)
    for j = 1:nrays
        dex = (temp(:,1) == j);
        for k = 1:numel(columnames)
            data.(columnames{k}){j,i} = temp(dex, k);
        end
    end
end

for k = 1:numel(columnames)
    tmp=data.(columnames{k});
    tmp =  stackWithPadND(tmp);
    tmp = reshape(tmp, [size(tmp,1), size(data.(columnames{k}) )]);     % L × nray × nbeams
    tmp=permute(tmp,[2 3 1]);
    data.(columnames{k}) =tmp;
end

columnames = {'reff_a','dPp_dV','dPt_dV','P_p','P_t','dP_dV','Pabs',...
    'jpar_p','jpar_t','Itor_p','Itor_t','effCD'};

[~,files] = get_files_sorted('Pabs_Icd_profiles*');

for i = 1:numel(files)
    S = importdata(files(i).name, ' ', 1);
    temp = S.data;
    for k = 1:numel(columnames)
        data.(columnames{k})(:,i) = temp(:, k);
    end
end
% for i=1:length(files)
%     temp = importdata(files(i).name,' ',1);
%     temp = temp.data;
%     data.reff(:,i) = temp(:,1);
%     data.dPpdV(:,i) = temp(:,2);
%     data.dPtdV(:,i) = temp(:,3);
%     data.Pp(:,i) = temp(:,4);
%     data.Pt(:,i) = temp(:,5);
%     data.dP0dV(:,i) = temp(:,6);
%     data.P0(:,i) = temp(:,7);
%     data.jcd_p(:,i) = temp(:,8);
%     data.jcd_t(:,i) = temp(:,9);
%     data.Itor_p(:,i) = temp(:,10);
%     data.Itor_t(:,i) = temp(:,11);
%     data.effCD(:,i) = temp(:,12);
% end
temp = importdata('nT_profiles',' ',2);
temp = temp.data;
data.profiles.reff = temp(:,1);
data.profiles.ne = temp(:,2);
data.profiles.te = temp(:,3);
data.profiles.Zeff = temp(:,4);
return;
end


function [idx,files] = get_files_sorted(expr)
% Get file list
listing = dir(expr);        % adjust extension or pattern as needed
names = {listing.name};
% Extract first integer in each filename
numTokens = regexp(names, '\d+', 'match');
% Convert to numeric
nums = cellfun(@(c) str2double(c{1}), numTokens);
% Sort by numeric value and reorder listing (and names if needed)
[~, idx] = sort(nums);
files = listing(idx);
end



% function A = stackWithPad(C)
%     n = numel(C);
%     L = cellfun(@numel, C);
%     M = max(L,[],'all');
%     A = NaN(M, n);
%     for i = 1:n
%             tmp=C{i};
%             out.(fn{k}) = cat(ndims(C)+1, vals{:});
%         A(1:L(i), i) = tmp;
%     end
%     A=reshape(A,M,)
% end

