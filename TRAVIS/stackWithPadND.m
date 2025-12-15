function A = stackWithPadND(C, padval)
% C: 1×N or N×1 cell array, each cell is a numeric array (any ndim)
% padval: optional pad value; default NaN for floating-point, 0 otherwise
if ~iscell(C) || isempty(C)
    error('C must be a non-empty cell array of numeric arrays.');
end
if ~all(cellfun(@isnumeric, C))
    error('All cells in C must be numeric.');
end

N = numel(C);

% Harmonize class (cast mixed classes to double)
classes = unique(cellfun(@class, C, 'UniformOutput', false));
if numel(classes) > 1
    C = cellfun(@double, C, 'UniformOutput', false);
    baseClass = 'double';
else
    baseClass = classes{1};
end

% Default pad value
if nargin < 2 || isempty(padval)
    if isfloat(feval(baseClass, 0))
        padval = NaN;
    else
        padval = 0;
    end
end

% Determine max dimensionality and target size across inputs
D = max(cellfun(@ndims, C),[],'all');             % next free dim will be D+1
sizes = zeros(N, D);
for i = 1:N
    si = size(C{i});
    sizes(i,1:numel(si)) = si;
    if numel(si) < D
        sizes(i,numel(si)+1:D) = 1;      % missing trailing dims are singleton
    end
end
target = max(sizes, [], 1);

% Preallocate output with padval
dimsOut = [target, N];
A = repmat(cast(padval, baseClass), dimsOut);

% Copy each input into its padded slot and stack along dim D+1
for i = 1:N
    idx = cell(1, D+1);
    for d = 1:D
        idx{d} = 1:sizes(i,d);
    end
    idx{D+1} = i;
    Ci = C{i};
    if ~strcmp(class(Ci), baseClass)
        Ci = cast(Ci, baseClass);
    end
    A(idx{:}) = Ci;
end
end
