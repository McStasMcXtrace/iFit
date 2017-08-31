function s = cat(varargin)
% concatenate structures

% first build up a cell array of structures
s_cell = {};

% handle input arguments
for index=1:numel(varargin)
  this = varargin{index};
  if isstruct(this) && numel(this) == 1, 
    % add a single structure
    s_cell{end+1} = this;
  elseif isstruct(this) && numel(this) > 1
    % add all items from a structure array
    for jj = 1:numel(this)
      s_cell{end+1} = this(jj);
    end
  end
end

% build up the field names and values
names = []; vals = {};
for index=1:numel(s_cell)
  names = [ names ; fieldnames(s_cell{index}) ];
  vals  = [ vals ; struct2cell(s_cell{index}) ];
end

vals  = vals(:);  % a column of values

% assemble the result

s = cell2struct(vals, names, 1);

%names = [fieldnames(struct1); fieldnames(struct2)];

% build up the structure values
%struct3 = cell2struct([struct2cell(struct1); struct2cell(struct2)], names, 1);

