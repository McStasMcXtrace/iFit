function outarray = iData(varargin)
% d = iData(a, ...) : iData class object constructor
%
% The iData objects are the containers where to import your data sets. It is then
% possible to Load, Plot, apply Math operators, Fit, and Save.
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
%   The syntax iData(iFunc object) evaluates the iFunc model using the iData
%     object axes, and returns the model value as an iData object.
%   The output argument is a single object or array of iData.
%   Input arguments may be more than one, or given as cells.
%
% When used with file names, compressed files may be used, as well as URLs and internal
%   anchor reference using the '#' character, such as in 'http://path/file.zip#Data'.
%
% Type <a href="matlab:doc(iData)">doc(iData)</a> to access the iFit/iData Documentation.
% Refer to <a href="matlab:helpwin iData/load">iData/load</a> for more information about loading data sets.
% iData is part of iFit http://ifit.mccode.org 
%
% examples:
%   d=iData('filename'); a=iData('http://filename.zip#Data');
%   d=iData(rand(10));
%
% Version: $Revision$
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
  a.Date         = clock;
  a.ModificationDate  = a.Date; % modification Date
  a.Command      = '';          % Matlab commands/history of the object
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
  a.Command      = { [ 'iData %% create ' a.Tag ] };
  a.Creator      = version(a);     % Creator (program) name
  outarray = [ outarray a ];
  return
else  % convert input argument into object
  if isa(varargin{1}, 'iData') && numel(varargin{1}) > 1
  % iData(iData)
    in = varargin{1};
    out = zeros(iData, numel(in), 1);
    parfor index=1:numel(in)
      out(index) = iData(in(index));        % check all elements
    end
    outarray = [ outarray out ];
    if nargout == 0 && ~isempty(inputname(1))
      assignin('caller',inputname(1),outarray)
    end
    return
  elseif ~isa(varargin{1}, 'iData') && isnumeric(varargin{1}) && length(varargin) > 1  % array -> iData
    % iData(x,y,..., signal)
    index = length(varargin);
    d = iData(varargin{index});  % last argument is the Signal
    
    % handle axes
    for k1=1:(index-1)
      % in plotting convention, X=2nd, Y=1st axis
      if     k1 <= 2 && ndims(d) >= 2, k2 = 3-k1; 
      else   k2 = k1; end
      set(d,    [ 'Data.Axis_' num2str(k1) ], varargin{k2});
      setaxis(d, k1, [ 'Axis_' num2str(k1) ], [ 'Data.Axis_' num2str(k1) ]);
      label(d, k1, inputname(k2));
    end
    % check in case the x,y axes have been reversed for dim>=2, then swap 1:2 axes in Signal
    if ndims(d)>=2 && isvector(getaxis(d, 1)) && isvector(getaxis(d, 2)) ...
                && length(getaxis(d, 1)) == size(get(d,'Signal'),2) ...
                && length(getaxis(d, 2)) == size(get(d,'Signal'),1) ...
                && length(getaxis(d, 1)) ~= length(getaxis(d, 2))
      s=get(d,'Signal'); set(d, 'Signal', s');
      disp([ 'iData: The Signal has been transposed to match the axes orientation in object ' d.Tag ' "' d.Title '".' ]);
    end
    if ~isempty(inputname(index))
        d.Label=[ inputname(index) ' (' class(varargin{index}) ')' ];
        label(d, 0, inputname(index));
        d.Title=inputname(index);
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
    elseif all(ishandle(in)) & numel(in)==1 % convert Handle Graphics Object
      % iData(figure handle)
      try 
        t = get(in,'DisplayName');
        if isempty(t), t=get(get(in,'Parent'),'DisplayName'); end
      catch
        t=[]; end
      if isempty(t), t=get(in,'Tag'); end
      if isempty(t), t=num2str(in); end
      if strcmp(get(in,'type'),'hggroup')
        t = [ 'figure ' t ];
        h = get(in,'Children');
        out = iData(h(1)); % fisrt item
        out = set(out,'Title', t);
        out = set(out, 'Label', t);
      elseif strcmp(get(in,'type'),'line')
        x = get(in,'xdata'); 
        y = get(in,'ydata'); 
        index = find(~isnan(x) & ~isnan(y));
        if length(index)~=numel(x), x = x(index); y=y(index); end
        c = get(in,'color');
        m = get(in,'marker');
        l = get(in,'linestyle');
        out=iData(x,y);
        try xl = get(get(in,'parent'),'XLabel'); xl=get(xl,'String'); catch 
            xl='x'; end; xlabel(out, xl);
        try yl = get(get(in,'parent'),'YLabel'); yl=[ get(yl,'String') ' ' ]; catch 
            yl=''; end;
        try tl = get(get(in,'parent'),'Title');  tl=[ get(tl,'String') ' ' ]; catch 
            tl=''; end;
        label(out,0,yl);
        t = [ 'line ' t ];
        out.Title = [ tl yl t ];
        out.DisplayName = t;
        out.Label=[ t ' marker ' m ' color ' num2str(c) ];
      elseif strcmp(get(in,'type'),'image')
        x = get(in,'xdata'); 
        y = get(in,'ydata');
        z = get(in,'cdata');
        t = [ 'image ' t ];
        out=iData(x,y,z);
        try xl = get(get(in,'parent'),'XLabel'); xl=get(xl,'String'); catch 
            xl='x'; end
        try yl = get(get(in,'parent'),'YLabel'); yl=get(yl,'String'); catch 
            yl='y'; end
        try zl = get(get(in,'parent'),'ZLabel'); zl=[ get(zl,'String') ' ' ]; catch 
            zl=''; end 
        try tl = get(get(in,'parent'),'Title');  tl=[ get(tl,'String') ' ' ]; catch 
            tl=''; end
        xlabel(out, xl); ylabel(out, yl); label(out, tl);
        out.Title = t;
        out.DisplayName = t;
        out.Label=t;
      elseif strcmp(get(in,'type'),'surface')
        x = get(in,'xdata'); 
        y = get(in,'ydata'); 
        z = get(in,'zdata'); 
        c = get(in,'cdata'); 
        % index=find(~isnan(x) & ~isnan(y) & ~isnan(z) & ~isnan(c)); 
        % if length(index)~=prod(size(x)), x = x(index); y=y(index); z=z(index); c=c(index); end
        l = get(in,'linestyle');
        if all(z == c)
          out=iData(x,y,z);
        else
          out=iData(x,y,z,c);
        end
        try xl = get(get(in,'parent'),'XLabel'); xl=get(xl,'String'); catch 
            xl='x'; end
        try yl = get(get(in,'parent'),'YLabel'); yl=get(yl,'String'); catch
            yl='y'; end
        try zl = get(get(in,'parent'),'ZLabel'); zl=[ get(zl,'String') ' ' ]; catch 
            zl=''; end 
        try tl = get(get(in,'parent'),'Title');  tl=[ get(tl,'String') ' ' ]; catch 
            tl=''; end
        xlabel(out, xl); ylabel(out, yl); label(out, tl);
        if all(z == c)
          t = [ tl zl t ];
        else
          if isempty(zl), zl='z'; end
          zlabel(out, zl);
          t = [ tl t ];
        end
        t = [ 'surface ' t ];
        out.Title = t;
        out.DisplayName = t;
        out.Label=[ t ' line ' l ];
      else
        h = [ findobj(in, 'type','line') ; findobj(in, 'type','surface') ; findobj(in, 'type','image')  ];
        out = [];
        for index=1:length(h)
          this_out = iData(h(index));
          if isempty(this_out.Title) && ~isempty(t)
            this_out.Title = t;
            this_out.Label = t;
            this_out.DisplayName = t;
          end
          if ~isempty(t), this_out.Source = t; end
          if  ~isscalar(get(this_out,'Signal'))
            out = [ out this_out ];
          end
        end
      end
    elseif isnumeric(in)
      % iData(x)
      out = iData_num2iData(in);    % convert scalar/vector/matrix to iData
    elseif iscell(in)
      % iData(cell)
      out = iData_cell2iData(in);   % convert cell/cellstr to cell(iData)
    elseif isa(in, 'iFunc')
      [signal, ax, name] = feval(in);
      if length(signal) == length(in.Parameters)
        [signal, ax, name] = feval(in, signal);
      end
      out = iData(ax{:}, signal);
      setalias(out, 'Error', 0);
    else
      iData_private_warning(mfilename, [ 'import of ' inputname(1) ' of class ' class(in) ' is not supported. Ignore.' ]);
      out = [];
    end
    if ~isempty(inputname(1)), in_name=[ inputname(1) ' ' ]; else in_name=''; end
    for index=1:numel(out)
      if numel(out) == 1 || ~isempty(out(index))
        if isempty(out(index).Source), out(index).Source = in_name; end
        if isempty(out(index).Title),  out(index).Title  = [ in_name ' (' class(in) ')' ]; end
        
        if isempty(out(index).Command)
          out(index) = iData_private_history(out(index), mfilename, in); 
        end
      end
    end
    out = iData_check(out); % private function
    if isa(in, 'iData') && nargout == 0 && ~isempty(inputname(1))
      assignin('caller',inputname(1),out);
    end
  end
end
outarray = [ outarray out ];

return

% ============================================================================
% iData_struct2iData: converts a structure into an iData
function b=iData_struct2iData(a)

  persistent fb

  if isempty(fb), fb=fieldnames(iData); end

  f  = fieldnames(a);
  b  = iData; 
  if isfield(a, 'Data')   % start by storing the raw Data
    b.Data = a.Data;
  end
  for index=1:length(f)
    if any(strcmp(f{index},fb))
      b = set(b,f{index}, a.(f{index}));
    end
  end
    
  if ~isfield(a, 'Data')   % store whole file content if possible.
    b.Data = a;
%  else
%    disp(['iData: warning: could not import all fields from structure.' ]);
  elseif isfield(a, 'Headers')
    b.Data.Headers = a.Headers;
    b=setalias(b, 'Headers', 'Data.Headers', 'Headers (text)' );
  end
  if isfield(a, 'Format')
    setalias(b, 'Format', a.Format);
  end
  if isfield(a, 'Command')
    b.Command = a.Command;
  end
  
  if isempty(b.Command), b.Command= cellstr('iData(<struct>)'); end
  
% ============================================================================
% iData_cell2iData: converts a cell into an iData array
function b=iData_cell2iData(a)
  b = [];
  for k=1:numel(a)
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
  setalias(b,'Signal','Data.Signal');
  if length(size(v)) > 2, v=v(:); end
  if numel(v) > 10, v=v(1:10); end
  v = mat2str(double(v)); 
  b.Command= cellstr([ 'iData(' v ')' ]);

% ============================================================================
function out = iData_check(in)
% make consistency checks on iData object

if numel(in) > 1
  out = zeros(iData, numel(in), 1);
  parfor index = 1:numel(in)
    out(index) = iData_check(in(index));
  end
  return
end

if iscell(in), in = in{1}; end

% update ModifDate
in.ModificationDate = clock;
% check type of fields
if iscellstr(in.Title)
  t = strcat(in.Title,';');
  in.Title=[ t{:} ];
end
if ~ischar(in.Title) 
  iData_private_warning(mfilename,['Title must be a char or cellstr in iData object ' in.Tag ' (' class(in.Title) '). Re-setting to empty.']);
  in.Title = '';
end
in.Title = strtrim(in.Title);
if ~ischar(in.Tag)
  iData_private_warning(mfilename,['Tag must be a char in iData object ' in.Tag ' "' in.Title '. Re-setting to a new Tad id.' ]);
  in = iData_private_newtag(in);
end
if ~ischar(in.Source)
  iData_private_warning(mfilename,['Source must be a char in iData object ' in.Tag ' "' in.Title '. Re-setting to pwd.' ]);
  in.Source = pwd;
end
if ~iscellstr(in.Command)
  in.Command = { in.Command };
end
if ~ischar(in.Date) && ~isnumeric(in.Date)
  iData_private_warning(mfilename,['Date must be a char/vector in iData object ' in.Tag ' "' in.Title '. Re-setting to "now".' ]);
  in.Date = clock;
end
if ~ischar(in.Creator)
  iData_private_warning(mfilename,['Creator must be a char in iData object ' in.Tag ' "' in.Title '. Re-setting to "iFit/iData".' ]);
  in.Creator = version(in);
end
if ~ischar(in.User)
  iData_private_warning(mfilename,['User must be a char in iData object ' in.Tag ' "' in.Title '. Re-setting to Matlab User.']);
  in.User = 'Matlab User';
end
% check if object.Data is numeric: make it a structure so that it is better organized
if isnumeric(in.Data) && ~isempty(in.Data)
  data = in.Data;
  in.Data = [];
  in.Data.Signal = data;
end

if ~isempty(in.Data) && isempty(getalias(in, 'Signal'))
  % get numeric fields sorted in descending size order
  [fields, types, dims] = findfield(in, '', 'numeric');
  if isempty(fields), 
    iData_private_warning(mfilename,['The iData object ' in.Tag ' "' in.Title '" contains no data at all ! (double/single/logical/int/uint)']);
  else
    fields_all = fields; dims_all=dims;
    % does this looks like a Signal ?
    
    if length(dims) > 1 % when similar sizes are encoutered, get the one which is not monotonic
      % list of 'biggest' fields
      maxdim=find(dims == dims(1)); maxdim2 = maxdim;
      % move 'error' and constant/monotonic down in the list
      for idx=1:length(maxdim)
        index=maxdim2(idx);
        x = get(in, fields{index});
        if ischar(x) || length(x) <= 1
          % this is a char: move to end of list
          maxdim([ end idx]) = maxdim([ idx end] );
          continue; 
        end
        x = diff(x(:));
        if all(x == x(1)) || all(x > 0) || ~isempty(strfind(lower(fields{index}), 'error'))
          % this is a constant/monotonic value or 'error'
          % move down in fields list
          maxdim([ end idx]) = maxdim([ idx end] );
        end
      end
      fields(maxdim2) = fields(maxdim);
    end
    
    % in case we have more than one choice, get the first one and error bars
    error_id = [];
    if length(dims) > 1 || iscell(fields)
      % do we have an 'error' which has same dimension ?
      for index=find(dims(:)' == dims(1))
        if index==1, continue; end % not the signal itself
        if ~isempty(strfind(lower(fields{index}), 'error'))
          error_id = fields{index};
        end
      end
      
      dims=dims(1);
      fields=fields{1};
    end

    % index: is the field and dimension index to assign the Signal
    if dims > 0
      disp([ 'iData: Setting Signal="' fields '" with length ' num2str(dims) ' in object ' in.Tag ' "' in.Title '".' ]);
      in = setalias(in,'Signal', fields);
      % search if there is a corresponding label (in Headers)
      lab = '';
      if isfield(in.Data, 'Headers')
        fields = fliplr(strtok(fliplr(fields), '.'));
        if isfield(in.Data.Headers, fields)
          in.Alias.Labels{1} = in.Data.Headers.(fields);
        end
      else
        label(in, 0, fields);
      end
      % assign potential 'error' bars
      if ~isempty(error_id)
        in = setalias(in,'Error', error_id);
      end
    end
    % look for vectors that may have the proper length as axes
    for index=1:ndims(in)
      if isempty(getaxis(in, num2str(index)))
        % search for a vector of length size(in, index)
        ax = find(dims_all == size(in, index));   % length of dim, or length(dim)+1
        if isempty(ax), ax = find(dims_all == size(in, index)+1); end
        ax = ax(~strcmp(getalias(in,'Signal'), fields_all(ax)));
        if length(ax) > 1; ax=ax(1); end
        if ~isempty(ax)
          val = get(in, fields_all{ax});
          if isvector(val) && ~strcmp(fields_all{ax},getalias(in,'Signal'))
            if length(val) == size(in, index) && min(val(:)) < max(val(:))
              in = setaxis(in, index, [ 'Axis_' num2str(index) ], fields_all{ax});
              found = 1;
            elseif length(val) == size(in, index)+1 && min(val(:)) < max(val(:))
              val = (val(1:(end-1)) + val(2:end))/2;
              in = setaxis(in, index, [ 'Axis_' num2str(index) ], val);
              found = 1;
            else found = 0;
            end
            if found == 1  % the axis could be found
              % search if there is a corresponding label (in Headers)
              if isfield(in.Data, 'Headers')
                fields=fliplr(strtok(fliplr(fields_all{ax}), '.'));
                if isfield(in.Data.Headers, fields)
                  in.Alias.Labels{index+1} = in.Data.Headers.(fields);
                end
              else
                label(in, index, fields_all{ax});
              end
              disp([ 'iData: Setting Axis{' num2str(index) '} ="' fields_all{ax} '" with length ' num2str(length(val)) ' in object ' in.Tag ' "' in.Title '".' ]);
            end
          end
          clear val
        else
          break; % all previous axes must be defined. If one misses, we end the search
        end
      end
    end 
  end
end

% check in case the x,y axes have been reversed for dim>=2, then swap 1:2 axes
if ndims(in)==2 && ~isempty(getaxis(in, '1')) && ~isempty(getaxis(in, '2')) ...
            && isvector(getaxis(in, 1)) && isvector(getaxis(in, 2)) ...
            && length(getaxis(in, 1)) == size(get(in,'Signal'),2) ...
            && length(getaxis(in, 2)) == size(get(in,'Signal'),1) ...
            && length(getaxis(in, 1)) ~= length(getaxis(in, 2))
  x1 = getaxis(in, '1');
  x2 = getaxis(in, '2');
  setaxis(in, 1, x2);
  setaxis(in, 2, x1);
  clear x1 x2
  disp([ 'iData: The object has been transposed to match the axes orientation in object ' in.Tag ' "' in.Title '".' ]);
end

% check aliases (valid ?) by calling setalias(in)
in = setalias(in);
% check axis (valid ?) by calling setaxis(in)
in = setaxis(in);

out = in;
