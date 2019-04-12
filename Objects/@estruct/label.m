function labl = label(this, varargin)
% b = label(s, alias, label) : Change iData label for a given alias/axis
%
%   @iData/label function to set/get labels
%     label(s, alias) returns the current label
%     label(s)        returns the object Label
%     label(s, rank, label)   sets the object label
%   The input iData object is updated if no output argument is specified.
%
% input:  s: object or array (iData)
%         alias: name of the alias or index of the axis (char/numeric)
%         label: new label (char/cellstr)
% output: b: object or array (iData)
% ex:     b=label(a,'x','new xlabel'); b=label(a,'x'); b=label(a, 1,'new xlabel');
%
% Version: $Date$ $Version$ $Author$
% See also iData, iData/plot, iData/xlabel, iData/ylabel, iData/zlabel, iDala/clabel

if numel(this) > 1
  if nargin < 3 % get labels (return them)
    labl=cell(size(this));
    for index=1:numel(this)
      labl{index} = label(this(index), varargin{:});
    end
  else          % set labels (update objects)
    for index=1:numel(this)
      label(this(index), varargin{:});
    end
    labl=this;
  end

  return
end

% handle different call syntax
if nargin == 1
  % label(a) -> get object Label property
  labl = this.Label;
  return
end

% get the alias we want to get/set (nargin == 2 or 3)
aliases = varargin{1};

% check input alias
if     ischar(aliases), aliases=cellstr(aliases); end
if ~iscellstr(aliases) && ~isnumeric(aliases)
  error([ mfilename ': Require a char/cellstr/numeric Alias name, but you gave me a ' class(aliases) ]);
end
if isnumeric(aliases), aliases=num2cell(aliases); end

% check input value (nargin==3)
if nargin >= 3
  values = varargin{2};
  if     ischar(values), values=cellstr(values); end
  if ~iscellstr(values)
    error([ mfilename ': Require a char/cellstr label, but you gave me a ' class(values) ]);
  end
else values = cell(size(aliases));
end

% we scan aliases
labl = {};
% a change of label should not trigger a check of Signal/axes.
if isfield(this.Private,'cache') && isfield(this.Private.cache,'check_requested')
  check = this.Private.cache.check_requested;
else
  check = false;
end
for index=1:numel(aliases)
  alias=aliases{index};

  l = label_single(this, alias, values{index}, nargin);
  if nargin == 2 % get
    if isempty(l)
      if strcmp(alias,'Signal')
        [l, alias] = label_single(this, '0', values{index}, nargin);
      elseif strcmp(alias,'0') || isequal(alias, 0)
        [l, alias] = label_single(this, 'Signal', values{index}, nargin);
      end
    end
    labl{end+1} = l;
  end

end
this.Private.cache.check_requested = check;
if nargin > 2, labl=this;
else
  if numel(labl) == 1, labl=labl{1}; end
end

% ------------------------------------------------------------------------------
function [labl, alias] = label_single(this, alias, value, n)
  % a single set/get for alias
  labl = '';
  if isnumeric(alias) || isfinite(str2double(alias)) % rank is given -> replace by corresponding alias
    tmp = getaxis(this, num2str(alias)); % this is the definition of axis rank
    if isempty(tmp) && alias == 0
      tmp = getaxis(this, 'Signal');
    end
  elseif isfield(this, alias)
    tmp = getaxis(this, alias);
  end
  if ischar(tmp) && isfield(this, tmp);
    alias = tmp;
  elseif isnumeric(alias)
    alias = sprintf('axis_%i', alias);
  end
  clear tmp
  if n == 2  % get
    if isfield(this, [ 'Labels.' alias ]) || isfield(this.Labels, alias) % is the alias label defined ?
         labl = subsref_single(this, [ 'Labels.' alias ]);  % ok we have it
    end                  % else invalid
  else            % set
    subsasgn_single(this, [ 'Labels.' alias ], validstr(value));
  end

% ------------------------------------------------------------------------------
function str=validstr(str)
  % validate a string as a single line
  if iscellstr(str), str=sprintf('%s;', str{:}); end
  if ~ischar(str), str=''; end % has been tested before (should not occur here)
  str=strrep(str(:)', sprintf('\n'), ';');
  index = find(str < 32 | str > 127);
  str(index) = ' ';
  str=strrep(str, '''', '''''');
  str=strrep(str,'\','/');
  if isempty(str), str=''; end
