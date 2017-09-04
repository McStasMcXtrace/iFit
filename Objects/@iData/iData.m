classdef iData
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
%       iData(1:100)
%     a string giving a file name to load. Use alternatively iData/load.
%       iData('filename_data')
%     a structure
%     a cell array which elements are imported separately
%     a iData object (updated if no output argument is specified).
%   The special syntax iData(x,y, .., c) creates an iData with
%     signal c and axes x,y, ...
%   The syntax iData(iFunc object, pars, axes) evaluates the iFunc model using the 
%     iData object axes, and returns the model value as an iData object.
%   The output argument is a single object or array of iData.
%   Input arguments may be more than one, or given as cells.
%
% When used with file names, compressed files may be used, as well as URLs and internal
%   anchor reference using the '#' character, such as in 'http://path/file.zip#Data'.
%
% Type <a href="matlab:doc(iData)">doc(iData)</a> to access the iFit/iData Documentation.
% Refer to <a href="matlab:doc iData/load">iData/load</a> for more information about loading data sets.
% iData is part of iFit http://ifit.mccode.org 
%
% examples:
%   d=iData('filename'); a=iData('http://filename.zip#Data');
%   d=iData(rand(10));
%
% Version: $Date$
% See also: iData, iData/load, methods, iData/setaxis, iData/setalias, iData/doc
% (c) E.Farhi, ILL. License: EUPL.

% object definition and converter
% EF 23/09/07 iData implementation
% ============================================================================
% internal functions
%   iData_struct2iData
%   iData_num2iData
%   iData_cell2iData
%   iData_check
%   load_clean_metadata

properties
  % create a new iData object
  Tag          = 0;           % unique ID
  Title        = '';          % Data title
  Source       = pwd;         % origin of data (filename/path)
  Creator      = [];          % Creator (program) name

  
  User         = '';          % User ID
  Date         = clock;
  ModificationDate  = ''; % modification Date
  Command      = '';          % Matlab commands/history of the object
  UserData     = '';          % user data storage area
  Label        = '';          % user label (color)
  DisplayName  = '';          % user name for handling data set as a variable
  
  % hidden fields
  Data         = [];          % Data storage area
  Alias        = {}
  
end % properties
  
methods
  function out = iData(varargin)
    % instantiate an iData object from input arguments
  
    user         = getenv('USER');
    if isempty(user), user = getenv('HOME'); end
    if isempty(user), user = getenv('TEMP'); end % gives User name on Windows systems
    
    out.Alias.Names  ={ ...
          'Signal','Error','Monitor'}; % (cell) Alias names
                                  % special alias: signal
                                  % special alias: errors on signal
                                  % special alias: monitor (statistical weight)
    out.Alias.Values ={'','',''};   % (cell) Alias values=string pointing to this Data/object or expr
    out.Alias.Labels ={...          % (cell) Alias labels/descriptions
          'Data Signal','Error on Signal','Monitor (weight)'};  
    out.Alias.Axis   ={};           % (cell) ordered list of axes names (Aliases)  
    
    out=iData_private_newtag(out);    % assign a new Tag/ID
    out.Command      = { [ 'iData %% create ' out.Tag ] };
    out.Creator      = version(out);    % Creator (program) name
    out.ModificationDate  = out.Date; % modification Date
    
    if nargin == 0, return; end

    if isa(varargin{1}, 'iData') && numel(varargin{1}) > 1
    % iData(iData array)
      out = varargin{1};
      for index=1:numel(out)
        out(index) = iData(out(index));        % check all elements
      end
      if nargout == 0 && ~isempty(inputname(1))
        assignin('caller',inputname(1),out)
      end
      return
    elseif ~isa(varargin{1}, 'iData') && isnumeric(varargin{1}) && length(varargin) > 1  % array -> iData
      % iData(x,y,..., signal)
      index = length(varargin);
      out   = iData(varargin{index});  % last argument is the Signal
      if ~isempty(inputname(index)) && any(~isfield(out, inputname(index)))
        out = setalias(out, inputname(index), 'Signal');
      end
      
      % handle axes
      for k1=1:(index-1)
        % in plotting convention, X=2nd, Y=1st axis
        % if     k1 <= 2 && ndims(out) >= 2, k2 = 3-k1; 
        % else   k2 = k1; end
        k2 = k1;
        out = set(out,    [ 'Data.Axis_' num2str(k1) ], varargin{k2});
        out = setaxis(out, k1, [ 'Axis_' num2str(k1) ], [ 'Data.Axis_' num2str(k1) ]);
        if ~isempty(inputname(k2))
          label(out, k1, inputname(k2));
          if ~isfield(out, inputname(k2))
            out = setalias(out, inputname(k2), [ 'Data.Axis_' num2str(k1) ]);
          end
        end
      end

      % check in case the x,y axes have been reversed for dim>=2, then swap 1:2 axes in Signal
      if ndims(out)>=2 && isvector(getaxis(out, 1)) && isvector(getaxis(out, 2)) ...
                  && length(getaxis(out, 1)) == size(get(out,'Signal'),2) ...
                  && length(getaxis(out, 2)) == size(get(out,'Signal'),1) ...
                  && length(getaxis(out, 1)) ~= length(getaxis(out, 2))
        s=get(out,'Signal'); out = set(out, 'Signal', s'); clear s
        disp([ 'iData: The Signal has been transposed to match the axes orientation in object ' out.Tag ' "' out.Title '".' ]);
      end
      
      if ~isempty(inputname(index))
          out.Label=[ inputname(index) ' (' class(varargin{index}) ')' ];
          out = label(out, 0, inputname(index));
          out.Title=inputname(index);
      end

      return
    elseif ischar(varargin{1}) % filename -> iData
    % iData('filename', ...)
      out = load(iData, varargin{:});        % load file(s) with additional arguments. Check included.
      return
    elseif isa(varargin{1}, 'iData')
      % iData(iData single)
      out = varargin{1};                     % just makes a check
    elseif isstruct(varargin{1})
      % iData(struct)
      out = iData_struct2iData(varargin{1}); % convert struct to iData
    elseif ~all(isempty(varargin{1})) && isnumeric(varargin{1}) && (numel(varargin{1}) > 1 || ~all(ishandle(varargin{1})))
      % iData(x)
      out = iData_num2iData(varargin{1});    % convert single scalar/vector/matrix to iData. No need for check.
      if ~isempty(inputname(1))
          out.Label=[ inputname(1) ' (' class(varargin{1}) ')' ];
          out = label(out, 0, inputname(1));
          out.Title=inputname(1);
      end
      return
    elseif ishandle(varargin{1}) % convert single Handle Graphics Object
      % iData(figure handle)
      out = iData_handle2iData(varargin{1});
    elseif iscell(varargin{1})
      % iData(cell)
      out = iData_cell2iData(varargin{1});   % convert cell/cellstr to cell(iData)
      return
    elseif isa(varargin{1}, 'iFunc')
      in = varargin{1};
      axes_in = varargin(3:end);
      out = [];
      
      for n_in = 1:numel(in)  % handle array of iFunc
        if numel(in) == 1, this_in = in;
        else               this_in = in(n_in); end
        
        [this_out, this_in] = iData_iFunc2iData(this_in, axes_in, varargin{2:end});
        
        if numel(in) == 1 || numel(in) == numel(this_in), in = this_in; 
        elseif numel(this_in) == 1
            in(n_in) = this_in; 
        end
        out = [ out this_out ];
      end % for n_in
      % update initial iFunc, if possible
      if ~isempty(inputname(1)) && numel(in) == numel(varargin{1})
         assignin('caller',inputname(1),in);
      end
      return
    else
      iData_private_warning(mfilename, [ 'import of ' inputname(1) ' of class ' class(varargin{1}) ' length ' mat2str(size(varargin{1})) ' is not supported. Ignore.' ]);
      out = [];
    end
      
    % check the object
    if ~isa(varargin{1}, 'iData')
      if ~isempty(inputname(1)), in_name=[ inputname(1) ' ' ]; else in_name=''; end
      for index=1:numel(out)
        if numel(out) == 1 || ~isempty(out(index))
          if isempty(out(index).Source), out(index).Source = in_name; end
          if isempty(out(index).Title),  out(index).Title  = [ in_name ' (' class(varargin{1}) ')' ]; end
          
          if isempty(out(index).Command)
            out(index) = iData_private_history(out(index), mfilename, varargin{1}); 
          end
        end
      end
    end
    
    out = iData_check(out); % private function
    
    if isa(varargin{1}, 'iData') && nargout == 0 && ~isempty(inputname(1))
      assignin('caller',inputname(1),out);
    end
      
  end % function iData

end % methods

end % classdef
  
% ============================================================================
% iData_cell2iData: converts a cell into an iData array
function b=iData_cell2iData(a)
  b = [];
  for k=1:numel(a)
    if numel(b) > 1, b=b(:); end
    c = iData(a{k});
    if numel(c) > 1, c=c(:); end
    b = [ b ; c ];
  end
  try
    b = reshape(b,size(a));
  end
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
end
% ============================================================================
% iData_num2iData: converts an iFunc model into an iData
function [out, this_in]=iData_iFunc2iData(this_in, axes_in, varargin)
  % evaluate the model 'this_in'
  % when some parameter are given as struct/cell and contain vectors, scans are done
  out = [];
  [signals, this_in, axs, names] = feval(this_in, varargin{:});
  
  if ~iscell(signals)
    signals = { signals };
    axs     = { axs };
    names   = { names };
  end
  % convert each of the returned data into an iData
  for i=1:numel(axs)
    signal = signals{i};
    ax     = axs{i};
    name   = names{i};
    if numel(this_in) > 1, this_model = this_in(i); else this_model=this_in; end
    
    if length(signal) == length(this_model.Parameters) % this was in fact a parameter guess...
      [signal, this_model, ax, name] = feval(this_model, signal, axes_in{:});
    end
    if isempty(signal), 
      iData_private_warning(mfilename, [ ': iFunc evaluation failed (empty value). Check axes and parameters.' ]);
      continue; 
    end
    
    % create iData object from signal and axes
    % check axes against signal size
    if numel(axes_in) == 1 && all(isnan(axes_in{1})), axes_in = []; end
    if isempty(axes_in) && ~isempty(ax), axes_in = ax; end
    for index=1:numel(axes_in)
      ax1 = axes_in{index};
      if index<numel(ax), ax2 = ax{index}; else ax2=[]; end
      if ~isempty(ax2)
        if iscell(ax1) && numel(ax1) == 1, ax1=ax1{1}; end
        if iscell(ax2) && numel(ax2) == 1, ax2=ax2{1}; end
        if isempty(ax1), ax1 = ax2;
        elseif ~isempty(find(isnan(ax1))), ax1 = ax2;
        elseif numel(ax2) == size(signal, index), ax1 = ax2;
        end
        axes_in{index} = ax1;
      end
    end
    this_out = iData(axes_in{:}, signal); % make it an iData

    % assign axes names
    if nargin > 2 % iData(iFunc,p,axes...)
      for index=1:numel(ax)
        if index+2 <= nargin && ~isempty(inputname(index+2))
          if numel(ax)>= 2 && index==1 && 0
            this_out=label(this_out,2,inputname(index+2)); 
          elseif numel(ax)>= 2 && index==2 && 0
            this_out=label(this_out,1,inputname(index+2));
          else
            this_out=label(this_out,index,inputname(index+2)); 
          end
        end
      end
    end
    this_out.Title = name;
    this_out.Label = name;
    this_out.DisplayName = name;
    setalias(this_out,'Error', 0);
    if ~isempty(this_model.ParameterValues)
        par_val = this_model.ParameterValues;
        pars    = this_model.Parameters;
      pars_out = cell2struct(num2cell(par_val(:)'), strtok(pars(:)'), 2);
      setalias(this_out,'Parameters', pars_out, [ name ' model parameters' ]);
    end
    setalias(this_out,'Model', this_model, this_model.Name);
    clear signal ax
    out = [ out this_out ];
  end % feval return arguments (can be a parameter scan)
  if numel(out) > 1 && numel(out) == numel(signals)
   out = reshape(out, size(signals));
  end
end
% ============================================================================

