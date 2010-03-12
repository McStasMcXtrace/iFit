function outarray = iData(varargin)
% d = iData(a, ...) : iData class object constructor
%
% Creates an iData object which contains data along with additional information.
%   An iData object may store any Data, and define its Signal, Error, Monitor, 
%     and Axes as aliases refering to e.g. parts of the Data.
%   The input argument 'a', when present, is converted into an iData object. It may be:
%     a scalar/vector/matrix
%     a string giving a file name to load. Use alternatively iData/load.
%     a structure
%     a cell array which elements are imported separately
%     a iData object (updated if no output argument is specified).
%   The special syntax iData(x,y, .., c) creates an iData with
%     signal c and axes x,y, ... where these are all numerics
%   The output argument is a single object or array of iData.
%   Input arguments may be more than one, or given as cells.
%
% examples:
%   d=iData('filename');
%   d=iData(rand(10));
%
% Version: $Revision: 1.14 $
% See also: iData, iData/load, methods

% object definition and converter
% EF 23/09/07 iData implementation
% ============================================================================
% internal functions
%   iData_struct2iData
%   iData_num2iData
%   iData_cell2iData
%   iData_check

outarray = [];

if nargin == 0  % create empty object
  % create a new iData object
  a.Title        = '';          % Data title
  
  a.Source       = pwd;         % origin of data (filename/path)
  a.Command      = cellstr('iData');          % Matlab commands/history of the object
  a.UserData     = '';          % user data storage area
  a.Label        = '';          % user label (color)
  a.DisplayName  = '';          % user name for handling data set as a variable
  a.Creator      = [];          % Creator (program) name
  user = getenv('USER');
  if isempty(user), user = getenv('HOME'); end
  if isempty(user), user = getenv('TEMP'); end % gives User name on Windows systems
  a.User         = user;        % User ID
  a.Data         =[];           % Data storage area
  a.Alias.Names  ={ ...
        'Signal','Error','Monitor'}; % (cell) Alias names
                                % special alias: signal
                                % special alias: errors on signal
                                % special alias: monitor (statistical weight)
  a.Alias.Values ={'','',''};   % (cell) Alias values=string pointing to this Data/object or expr
  a.Alias.Labels ={...          % (cell) Alias labels/descriptions
        'Data Signal','Error on Signal','Monitor (weight)'};  
  a.Alias.Axis   ={};           % (cell) ordered list of axes names (Aliases)  
  
  a=iData_private_newtag(a);     
  a.ModificationDate  = a.Date; % modification Date

  % create the object
  a = class(a, 'iData');
  a.Command      = cellstr([ 'iData %% create ' a.Tag ]);
  a.Creator      =  version(a);     % Creator (program) name
  outarray = [ outarray a ];
  return
else  % convert input argument into object
  if isa(varargin{1}, 'iData') & length(varargin{1}) > 1
    in = varargin{1};
    for index=1:length(in)
      out(index) = iData(in(index));        % check all elements
    end
    outarray = [ outarray out ];
    if nargout == 0 & length(inputname(1))
      assignin('caller',inputname(1),outarray)
    end
    return
  elseif isnumeric(varargin{1}) & length(varargin) > 1  % array -> iData
    index=1;
    while index <= length(varargin)
      if isnumeric(varargin{index})
        for j1=index:length(varargin)
          if ~isnumeric(varargin{j1})
            j1 = max(index, j1-1); 
            break; 
          end
        end
        % numerics are from i to j1
        d = iData(varargin{j1});
        for k1=1:(j1-index)
          d = set(d, [ 'Data.axis' num2str(k1) ], varargin{index+k1-1});
          d = setaxis(d, k1, [ 'axis' num2str(k1) ], [ 'Data.axis' num2str(k1) ]);
        end
        index = j1;
        outarray = [ outarray d ];
      else
        outarray = [ outarray iData(varargin{i}) ];
      end
      index = index+1;
    end % while
    return
  elseif ischar(varargin{1}) & length(varargin) > 1 % filename -> iData
    out = load(iData, varargin{:});        % load file(s) with additional arguments

  else
    in = varargin{1};
    if ischar(in)
      out = load(iData, in);        % load file(s)
    elseif isa(in, 'iData')
      out = in;                     % just makes a check
    elseif isstruct(in)
      out = iData_struct2iData(in); % convert struct to iData
    elseif ishandle(in)             % convert Handle Graphics Object
      if strcmp(get(in,'type'),'hggroup')
        t = get(in,'DisplayName');
        if isempty(t), t=get(in,'Tag'); end
        h = get(in,'Children');
        out = iData(h(1)); % fisrt item
        out = set(out,'Title', t);
        out = set(out, 'DisplayName', t);
      elseif strcmp(get(in,'type'),'line')
        x = get(in,'xdata'); x = x(find(~isnan(x)));
        y = get(in,'ydata'); y = y(find(~isnan(y)));
        t = get(in,'DisplayName');
        if isempty(t), t=get(in,'Tag'); end
        c = get(in,'color');
        m = get(in,'marker');
        l = get(in,'linestyle');
        out=iData(x,y);
        out.Title = t;
        out.DisplayName = t;
        out.Label=[ 'line ' l ' marker ' m ' color ' num2str(c) ];
      elseif strcmp(get(in,'type'),'surface')
        x = get(in,'xdata'); x = x(find(~isnan(x)));
        y = get(in,'ydata'); y = y(find(~isnan(y)));
        z = get(in,'zdata'); z = z(find(~isnan(z)));
        c = get(in,'cdata'); c = c(find(~isnan(c)));
        l = get(in,'linestyle');
        t = get(in,'DisplayName');
        if isempty(t), t=get(in,'Tag'); end
        if all(z == c)
          out=iData(x,y,z);
        else
          out=iData(x,y,z,c);
        end
        out.Title = t;
        out.DisplayName = t;
        out.Label=[ 'surface ' t ' line ' l ];
      else
        h = [ findobj(in, 'type','line') findobj(in, 'type','surface') ];
        try t = get(in,'DisplayName'); catch t=[]; end
        if isempty(t)
          try t = get(in,'Tag'); catch t=[]; end
        end
        out = [];
        for index=1:length(h)
          this_out = iData(h(index));
          if isempty(this_out.Title) && ~isempty(t)
            this_out.Title = t;
            this_out.DisplayName = t;
          end
          if ~isempty(t), this_out.Source = t; end
          out = [ out this_out ];
        end
      end
    elseif isnumeric(in)
      out = iData_num2iData(in);    % convert scalar/vector/matrix to iData
    elseif iscell(in)
      out = iData_cell2iData(in);   % convert cell/cellstr to cell(iData)
    else
      iData_private_warning(mfilename, [ 'import of ' inputname(1) ' of class ' class(in) ' is not supported. Ignore.' ]);
      out = [];
    end
    if length(inputname(1)), inmame=inputname(1); else inmame=''; end
    for index=1:length(out)
      if length(out) == 1 | ~isempty(out(index))
        if isempty(out(index).Source), out(index).Source = inmame; end
        if isempty(out(index).Title),  out(index).Title  = [ inmame ' ' class(in) ' import into iData ' ]; end
        
        if isempty(out(index).Command)
        	out(index) = iData_private_history(out(index), mfilename, in); 
        end
        out(index) = iData_check(out(index));  % private function
      end
    end
    if isa(in, 'iData') & nargout == 0
      assignin('caller',inputname(1),out);
    end
  end
end
outarray = [ outarray out ];

return

% ============================================================================
% iData_struct2iData: converts a structure into an iData
function b=iData_struct2iData(a)

  f  = fieldnames(a);
  b  = iData; 
  fb = fieldnames(b);
  if isfield(a, 'Data')   % start by storing the raw Data
    b.Data = a.Data;
  end
  for index=1:length(f)
    if ~isempty(strmatch(f{index},fb,'exact'))
      b = set(b,f{index}, getfield(a,f{index}));
    end
  end
    
  if ~isfield(a, 'Data')   % store whole file content if possible.
    b.Data = a;
%  else
%    disp(['iData: warning: could not import all fields from structure.' ]);

  if isempty(b.Command), b.Command= cellstr('iData(<struct>)'); end
  end
  
% ============================================================================
% iData_cell2iData: converts a cell into an iData cell
function b=iData_cell2iData(a)
  b = [];
  for k=1:length(a(:))
    b = [ b iData(a{k}) ];
  end
  b = reshape(b,size(a));

% ============================================================================
% iData_num2iData: converts a numeric into an iData
function b=iData_num2iData(v)
  b=iData;
  b.Data.Signal = v;
  if length(size(v)) > 2, v=v(:); end
  if numel(v) > 100, v=v(1:100); end
  v = mat2str(double(v)); 
  b.Command= cellstr([ 'iData(' v ')' ]);

% ============================================================================
function out = iData_check(in)
% make consistency checks on iData object

if length(in) > 1
  for index = 1:length(in)
    out(index) = iData_check(in(index));
  end
  return
end

if iscell(in), in = in{1}; end

% update ModifDate
in.ModificationDate = datestr(now);
% check type of fields
if ~ischar(in.Title) & ~iscellstr(in.Title)
  iData_private_warning(mfilename,'Title must be a char or cellstr');
  in.Title = '';
end
if ~ischar(in.Tag)
  iData_private_warning(mfilename,'Tag must be a char');
  in = iData_private_newtag(in);
end
if ~ischar(in.Source)
  iData_private_warning(mfilename,'Source must be a char');
  in.Source = '';
end
if ~ischar(in.Command) & ~iscellstr(in.Command)
  iData_private_warning(mfilename,'Command must be a char or cellstr');
  in.Command = cellstr('');
endif ~ischar(in.Date)
  iData_private_warning(mfilename,'Date must be a char');
  in.Date = datestr(now);
end
if ~ischar(in.Creator)
  iData_private_warning(mfilename,'Creator must be a char');
  in.Creator = version(in);
endif ~ischar(in.User)
  iData_private_warning(mfilename,'User must be a char');
  in.User = 'Matlab User';
end
if isempty(in.Data)
  in = setalias(in, getalias(in));
% if signal is invalid, set signal to biggest field link
elseif isempty(getalias(in, 'Signal'))
  [fields, types, dims] = findfield(in);
  index=strmatch('double', types, 'exact');
  if isempty(index), index=strmatch('single', types, 'exact'); end
  if isempty(index), index=strmatch('logical', types, 'exact'); end
  if isempty(index), index=strmatch('uint', types); end
  if isempty(index), 
    iData_private_warning(mfilename,['The iData object ' in.Tag ' contains no data at all ! (double/single/logical)']);
  else
    fields = fields(index); % get all field names containing double data
    dims = dims(index);
    [dummy, index] = sort(dims);
    index=index(end);
    if dummy(index) > 0
      disp([ 'iData: Setting the Signal of ' in.Tag ' to the biggest numerical field ' fields{index} ' with length ' num2str(dummy(end)) '.' ]);
      in = setalias(in,'Signal', fields{index});
    end
  end
end
% check aliases (valid ?) by calling setalias(in)
in = setalias(in);
% check axis (valid ?) by calling setaxis(in)
in = setaxis(in);

% and make it an iData object
out = class(struct(in), 'iData');


