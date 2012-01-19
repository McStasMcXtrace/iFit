function outarray = iData(varargin)
% d = iData(a, ...) : iData class object constructor
%
% Creates an iData object which contains data along with additional information.
%   An iData object may store any Data, and define its Signal, Error, Monitor, 
%     and Axes as aliases refering to e.g. parts of the Data.
%   The input argument 'a', is converted into an iData object. It may be:
%     a scalar/vector/matrix
%     a string giving a file name to load. Use alternatively iData/load.
%     a structure
%     a cell array which elements are imported separately
%     a iData object (updated if no output argument is specified).
%   The special syntax iData(x,y, .., c) creates an iData with
%     signal c and axes x,y, ... where these are all numerics with 'x'
%     for the columns (2nd axis rank), 'y' for the rows (1st axis rank), ...
%   The output argument is a single object or array of iData.
%   Input arguments may be more than one, or given as cells.
%
% When used with file names, compressed files may be used, as well as URLs and internal
%   anchor reference using the '#' character, such as in 'http://path/file.zip#Data'.
% Refer to the iData/load command for more information.
%
% examples:
%   d=iData('filename'); a=iData('http://filename.zip#Data');
%   d=iData(rand(10));
%
% Version: $Revision: 1.34 $
% See also: iData, iData/load, methods, iData/setaxis, iData/setalias, iData/doc

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
  a.Tag          = 0;           % unique ID
  a.Title        = '';          % Data title
  a.Source       = pwd;         % origin of data (filename/path)
  a.Creator      = [];          % Creator (program) name
  user = getenv('USER');
  if isempty(user), user = getenv('HOME'); end
  if isempty(user), user = getenv('TEMP'); end % gives User name on Windows systems
  a.User         = user;        % User ID
  a.Date         = datestr(now);;
  a.ModificationDate  = a.Date; % modification Date
  a.Command      = cellstr('iData');          % Matlab commands/history of the object
  a.UserData     = '';          % user data storage area
  a.Label        = '';          % user label (color)
  a.DisplayName  = '';          % user name for handling data set as a variable
  
  % hidden fields
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

  % create the object
  a=iData_private_newtag(a); 
  a = class(a, 'iData');
  a.Command      = cellstr([ 'iData %% create ' a.Tag ]);
  a.Creator      = version(a);     % Creator (program) name
  outarray = [ outarray a ];
  return
else  % convert input argument into object
  if isa(varargin{1}, 'iData') & length(varargin{1}) > 1
  % iData(iData)
    in = varargin{1};
    for index=1:length(in)
      out(index) = iData(in(index));        % check all elements
    end
    outarray = [ outarray out ];
    if nargout == 0 & length(inputname(1))
      assignin('caller',inputname(1),outarray)
    end
    return
  elseif ~isa(varargin{1}, 'iData') && isnumeric(varargin{1}) & length(varargin) > 1  % array -> iData
    % iData(x,y,..., signal)
    index = length(varargin);
    d = iData(varargin{index});  % last argument is the Signal
    % handle axes
    for k1=1:(index-1)
      % in plotting convention, X=2nd, Y=1st axis
      if     k1 <= 2 & ndims(d) >= 2, k2 = 3-k1; 
      else   k2 = k1; end
      set(d,    [ 'Data.Axis_' num2str(k1) ], varargin{k2});
      setaxis(d, k1, [ 'Axis_' num2str(k1) ], [ 'Data.Axis_' num2str(k1) ]);
    end
    % check in case the x,y axes have been reversed for dim>=2, then swap 1:2 in Signal
    if ndims(d)>=2 && isvector(getaxis(d, 1)) && isvector(getaxis(d, 2)) ...
                && length(getaxis(d, 1)) == size(get(d,'Signal'),2) ...
                && length(getaxis(d, 2)) == size(get(d,'Signal'),1)
      s=get(d,'Signal'); set(d, 'Signal', s');
      disp([ 'iData: The Signal has been transposed to match the axes orientation in object ' d.Tag ' "' d.Title '".' ]);
    end
      
    outarray = [ outarray d ];
    return
  elseif ischar(varargin{1}) & length(varargin) > 1 % filename -> iData
  % iData('filename', ...)
    out = load(iData, varargin{:});        % load file(s) with additional arguments
  else
    in = varargin{1};
    if ischar(in)
      % iData('filename')
      out = load(iData, in);        % load file(s)
      if ~isa(out, 'iData'), outarray=out; return; end
    elseif isa(in, 'iData')
      % iData(iData)
      out = in;                     % just makes a check
    elseif isstruct(in)
      % iData(struct)
      out = iData_struct2iData(in); % convert struct to iData
    elseif ishandle(in)             % convert Handle Graphics Object
      % iData(figure handle)
      if strcmp(get(in,'type'),'hggroup')
        try t = get(in,'DisplayName'); catch t=[]; end
        if isempty(t), t=get(in,'Tag'); end
        h = get(in,'Children');
        out = iData(h(1)); % fisrt item
        out = set(out,'Title', t);
        out = set(out, 'DisplayName', t);
      elseif strcmp(get(in,'type'),'line')
        x = get(in,'xdata'); 
        y = get(in,'ydata'); 
        index = find(~isnan(x) & ~isnan(y));
        if length(index)~=prod(size(x)), x = x(index); y=y(index); end
        try t = get(in,'DisplayName'); catch t=[]; end
        if isempty(t), t=get(in,'Tag'); end
        c = get(in,'color');
        m = get(in,'marker');
        l = get(in,'linestyle');
        out=iData(x,y);
        try xl = get(get(in,'parent'),'XLabel'); xl=get(xl,'String'); catch xl='x'; end; xlabel(out, xl);
        try yl = get(get(in,'parent'),'YLabel'); yl=[ get(yl,'String') ' ' ]; catch yl=''; end;
        try tl = get(get(in,'parent'),'Title');  tl=[ get(tl,'String') ' ' ]; catch tl=''; end;
        label(out,0,yl);
        out.Title = [ tl yl t ];
        out.DisplayName = t;
        out.Label=[ 'line ' l ' marker ' m ' color ' num2str(c) ];
      elseif strcmp(get(in,'type'),'image')
        x = get(in,'xdata'); 
        y = get(in,'ydata');
        z = get(in,'cdata');
        t = get(in,'Tag');
        out=iData(x,y,z);
        try xl = get(get(in,'parent'),'XLabel'); xl=get(xl,'String'); catch xl='x'; end
        try yl = get(get(in,'parent'),'YLabel'); yl=get(yl,'String'); catch yl='y'; end
        try zl = get(get(in,'parent'),'ZLabel'); zl=[ get(zl,'String') ' ' ]; catch zl=''; end 
        try tl = get(get(in,'parent'),'Title');  tl=[ get(tl,'String') ' ' ]; catch tl=''; end
        xlabel(out, xl); ylabel(out, yl); label(out, tl);
        out.Title = t;
        out.DisplayName = t;
        out.Label=[ 'image ' t ];
      elseif strcmp(get(in,'type'),'surface')
        x = get(in,'xdata'); 
        y = get(in,'ydata'); 
        z = get(in,'zdata'); 
        c = get(in,'cdata'); 
        % index=find(~isnan(x) & ~isnan(y) & ~isnan(z) & ~isnan(c)); 
        % if length(index)~=prod(size(x)), x = x(index); y=y(index); z=z(index); c=c(index); end
        l = get(in,'linestyle');
        t = get(in,'DisplayName');
        if isempty(t), t=get(in,'Tag'); end
        if all(z == c)
          out=iData(x,y,z);
        else
          out=iData(x,y,z,c);
        end
        try xl = get(get(in,'parent'),'XLabel'); xl=get(xl,'String'); catch xl='x'; end
        try yl = get(get(in,'parent'),'YLabel'); yl=get(yl,'String'); catch yl='y'; end
        try zl = get(get(in,'parent'),'ZLabel'); zl=[ get(zl,'String') ' ' ]; catch zl=''; end 
        try tl = get(get(in,'parent'),'Title');  tl=[ get(tl,'String') ' ' ]; catch tl=''; end
        xlabel(out, xl); ylabel(out, yl); label(out, tl);
        if all(z == c)
          t = [ tl zl t ];
        else
          if isempty(zl), zl='z'; end
          zlabel(out, zl);
          t = [ tl t ];
        end
        out.Title = t;
        out.DisplayName = t;
        out.Label=[ 'surface ' t ' line ' l ];
      else
        h = [ findobj(in, 'type','line') findobj(in, 'type','surface') findobj(in, 'type','image')  ];
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
      % iData(x)
      out = iData_num2iData(in);    % convert scalar/vector/matrix to iData
    elseif iscell(in)
      % iData(cell)
      out = iData_cell2iData(in);   % convert cell/cellstr to cell(iData)
    else
      iData_private_warning(mfilename, [ 'import of ' inputname(1) ' of class ' class(in) ' is not supported. Ignore.' ]);
      out = [];
    end
    if length(inputname(1)), inmame=[ inputname(1) ' ' ]; else inmame=''; end
    for index=1:length(out)
      if length(out) == 1 | ~isempty(out(index))
        if isempty(out(index).Source), out(index).Source = inmame; end
        if isempty(out(index).Title),  out(index).Title  = [ inmame class(in) ]; end
        
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
    if any(strcmp(f{index},fb))
      b = set(b,f{index}, getfield(a,f{index}));
    end
  end
    
  if ~isfield(a, 'Data')   % store whole file content if possible.
    b.Data = a;
%  else
%    disp(['iData: warning: could not import all fields from structure.' ]);
  end
  % case of importation from an iFunc model structure
  if isfield(a, 'Type') && strcmp(a.Type,'iFit fitting function')
    if ~isempty(a.Guess)
      try
      b.Data.Parameters=cell2struct(num2cell(a.Guess), a.Parameters, 2);
      setalias(b, 'Parameters', 'Data.Parameters',[ a.function ' model parameters' ]);
      end
    end
    setalias(b, 'Signal', 'Data.Values');
    set(b, 'Error', 0);
    title(b, a.function);
    for index=1:length(a.Axes)
      setaxis(b, index, [ 'Data.Axes{' num2str(index) '}' ]);
      b.Title=a.Name;
    end
  end
  if isempty(b.Command), b.Command= cellstr('iData(<struct>)'); end
  
% ============================================================================
% iData_cell2iData: converts a cell into an iData cell
function b=iData_cell2iData(a)
  b = [];
  for k=1:length(a(:))
    b = [ b iData(a{k}) ];
  end
  try
      b = reshape(b,size(a));
  end

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
if iscellstr(in.Title)
  t = strcat(in.Title,';');
  in.Title=[ t{:} ];
end
if ~ischar(in.Title) 
  iData_private_warning(mfilename,['Title must be a char or cellstr in iData object ' in.Tag ' "' in.Title ]);
  in.Title = '';
end
if ~ischar(in.Tag)
  iData_private_warning(mfilename,['Tag must be a char in iData object ' in.Tag ' "' in.Title ]);
  in = iData_private_newtag(in);
end
if ~ischar(in.Source)
  iData_private_warning(mfilename,['Source must be a char in iData object ' in.Tag ' "' in.Title ]);
  in.Source = '';
end
if ~ischar(in.Command) & ~iscellstr(in.Command)
  iData_private_warning(mfilename,['Command must be a char or cellstr in iData object ' in.Tag ' "' in.Title ]);
  in.Command = cellstr('');
end
if ~ischar(in.Date)
  iData_private_warning(mfilename,['Date must be a char in iData object ' in.Tag ' "' in.Title ]);
  in.Date = datestr(now);
end
if ~ischar(in.Creator)
  iData_private_warning(mfilename,['Creator must be a char in iData object ' in.Tag ' "' in.Title ]);
  in.Creator = version(in);
end
if ~ischar(in.User)
  iData_private_warning(mfilename,['User must be a char in iData object ' in.Tag ' "' in.Title ]);
  in.User = 'Matlab User';
end
% check if object.Data is numeric: make it a structure so that it is better organized
if isnumeric(in.Data) && ~isempty(in.Data)
  data = in.Data;
  in.Data = [];
  in.Data.Signal = data;
end

if isempty(in.Data)
  in = setalias(in, getalias(in));
% if signal is invalid, set signal to biggest field link
elseif isempty(getalias(in, 'Signal'))
  [fields, types, dims] = findfield(in);
  index=find(strcmp('double', types));
  index=[ index ; find(strcmp('single',  types)) ];
  index=[ index ; find(strcmp('logical', types)) ];
  index=[ index ; find(strncmp('uint', types, 4)) ];
  if isempty(index), 
    iData_private_warning(mfilename,['The iData object ' in.Tag ' "' in.Title '" contains no data at all ! (double/single/logical)']);
  else
    fields = fields(index); % get all field names containing double data
    dims = dims(index);
    % now we get the largest dim, and reduce the search list
    index = find(dims == max(dims));
    dims  = dims(index);
    fields= fields(index);
    % does this look like a Signal ?
    if length(dims) > 1 % when similar sizes are encoutered, get the one which is not monotonic
      for index=1:length(dims)
        x = get(in, fields{index});
        if ischar(x) || length(x) <= 1, continue; end
        x = diff(x(:));
        if all(x == x(1)) || all(x > 0) % this is a constant/monotonic step axis
          continue;
        elseif any(x < 0) && any(x> 0)
          dims  = dims(index);
          fields= fields{index};
          break;
        end
      end
    end   
    % in case we still have more than one choice, get the last one
    if length(dims) > 1 || iscell(fields)
      dims=dims(end);
      fields=fields{end};
    end
    % index: is the field and dimension index to assign the Signal
    if dims > 0
      disp([ 'iData: Setting Signal="' fields '" with length ' num2str(dims) ' in object ' in.Tag ' "' in.Title '".' ]);
      in = setalias(in,'Signal', fields);
    end
  end
end
% check in case the x,y axes have been reversed for dim>=2, then swap 1:2 in Signal
if ndims(in)==2 && ~isempty(getaxis(in, '1')) && ~isempty(getaxis(in, '2')) ...
            && isvector(getaxis(in, 1)) && isvector(getaxis(in, 2)) ...
            && length(getaxis(in, 1)) == size(get(in,'Signal'),2) ...
            && length(getaxis(in, 2)) == size(get(in,'Signal'),1)
  s=get(in,'Signal'); set(in, 'Signal', s');
  disp([ 'iData: The Signal has been transposed to match the axes orientation in object ' in.Tag ' "' in.Title '".' ]);
end
    
% check aliases (valid ?) by calling setalias(in)
in = setalias(in);
% check axis (valid ?) by calling setaxis(in)
in = setaxis(in);

% and make it an iData object
out = class(struct(in), 'iData');

