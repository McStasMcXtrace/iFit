function labl = label(this, rank, lab)
% b = label(s, alias, label) : Change iData label for a given alias/axis
%
%   @iData/label function to set/get labels
%     label(s, alias) returns the current label
%   The input iData object is updated if no output argument is specified.
%
% input:  s: object or array (iData)
%         alias: name of the alias or index of the axis (char/numeric)
%         label: new label (char/cellstr)
% output: b: object or array (iData)
% ex:     b=label(a,'x','new xlabel'); b=label(a,'x'); b=label(a, 1,'new xlabel');
%
% Version: $Revision: 1.10 $
% See also iData, iData/plot, iData/xlabel, iData/ylabel, iData/zlabel, iDala/clabel

if nargin < 2, rank=[]; end
if nargin < 3, lab=[]; end

if numel(this) > 1
  if nargin < 3
    labl=cell(size(this));
    for index=1:numel(this)
      labl{index} = label(this(index), rank, lab);
    end
  else
    for index=1:numel(this)
      this(index) = label(this(index), rank, lab);
    end
    labl=this;
    if nargout == 0 && ~isempty(inputname(1))
      assignin('caller',inputname(1),this);
    end
  end
  
  return
end

% label(a)
if isempty(rank), labl=this.Label; return;
else              labl = '';
end

alias = '';
if iscellstr(rank)
  rank = rank{1};
end
if ischar(rank)  % label(a, '1', ...)
  if ~isnan(str2double(rank))  % e.g. alias = '1'
    rank = str2double(rank);  % now numeric axis index
  else            % label(a, 'X', ...)
    % alias is an Alias reference, search for it in the Axes
    if strcmp(rank, 'Signal'), rank = 0;
    elseif any(strcmp(rank, this.Alias.Axis)), rank  = find(strcmp(rank, this.Alias.Axis)); % now numeric axis index or empty
    elseif any(strcmp(rank, this.Alias.Axis)), alias = find(strcmp(rank, this.Alias.Names)); rank = -1;  % this is an alias already, NOT an axis
    else alias = rank; rank = -1; % new alias to create
    end
  end
end

if isempty(rank),    return; end  % axis alias not found
%if ~isnumeric(rank), return; end

% label(a, 1, ...): get the definition of the axis (alias name)
value = '';
if ~isempty(rank) && rank == 0
  alias = 'Signal';
elseif ~isempty(rank) && rank > 0 && rank <= length(this.Alias.Axis)
  alias = this.Alias.Axis{rank};  % name of the alias/axis (or direct value storage)
  if ~ischar(alias), value = alias; alias = ''; end
end

% does the axis alias already exists ? if not create it so that we can set its label
if nargin == 3 && ~isempty(lab)
  if rank > 0
    % create an axis
    if isempty(value), value = getaxis(this, rank); end % probably the default axis
    setaxis(this, rank, value);
    alias = getaxis(this, num2str(rank));
  else
    try
    value = get(this, alias); 
    catch
    end
    if isempty(value), value = 0; end
    setalias(this, alias, value, lab);
  end
end

% searches the axis definition in the alias names
alias_num   = find(strcmpi(alias, this.Alias.Names));
if ~isempty(alias_num)
  labl = this.Alias.Labels{alias_num};
end

if nargin == 2 || isempty(lab)  % return the current label value
  return;
end

% now we have nargin == 3
if ~isempty(alias_num)  % the alias already exists: store the new alias value
  lab = regexprep(lab,'\s+',' '); % remove duplicated spaces
  this.Alias.Labels{alias_num} = lab;
  if ~isempty(value)
    this.Alias.Values{alias_num} = value;
  end
end

this = iData_private_history(this, mfilename, this, rank, lab);
labl = this;

if nargout == 0 && ~isempty(inputname(1))
  assignin('caller',inputname(1),this);
end

