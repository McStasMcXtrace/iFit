function s = setaxis(s,varargin)
% SETAXIS set axis definition or value in object.
%   SETAXIS(a, rank, value) Set axis value (and follow aliases). This is equivalent to
%   the syntax a{rank}=value. The axis rank 0 corresponds with the Signal/Monitor value.
%   The axis of rank 1 corresponds with rows, 2 with columns, 3 with pages, etc.
%
%   SETAXIS(a, 'Signal',value) Set the Signal/Monitor value, equivalent to
%   a{0}=value and setaxis(a, 0, value).
%
%   SETAXIS(a, 'Error', value) Set the Error/Monitor value.
%
%   SETAXIS(a, 'rank', 'alias') Set axis definition (alias). The axis of rank 0
%   corresponds with the Signal definition. The 'value' can be specified as
%  'biggest' to indicate the numeric biggest array, as in findfield.
%
%   An alias is a string/char which allows to link to internal or external links
%   as well as evaluated expression, with the following syntax cases:
%     'field'                           a simple link to an other property 'field'
%     'field1.field2...'                a nested link to an other property
%     'file://some_file_path'           a local file URL
%     'http://some_distant_resource'    an HTTP URL (proxy settings may need to be set)
%     'https://some_distant_resource'   an HTTPS URL (proxy settings may need to be set)
%     'ftp://some_distant_resource'     an FTP URL (proxy settings may need to be set)
%     'matlab: some_expression'         some code to evaluate. 'this' refers to the object itself
%
%   File and URL can refer to compressed resources (zip, gz, tar, Z) which are
%   extracted on-the-fly. In case the URL/file resource contains 'sections', a
%   search token can be specified with syntax such as 'file://filename#token'.
%
%   SETAXIS(a) check axes, Signal, Monitor, Error.
%
% Example: s=estruct('x',1:10,'y',1:20, 'data',rand(10,20)); setaxis(s,1,'x'); isnumeric(getaxis(s,1))
% Version: $Date$ $Version$ $Author$
% See also estruct, fieldnames, findfield, isfield, set, get, getalias, setalias,
%   getaxis, setaxis

  if nargin == 1, s = axescheck(s); return; end

  % handle array of struct
  if numel(s) > 1
    for index=1:numel(s)
      s(index) = setaxis(s(index), varargin{:});
    end
    return
  end

  m = []; % will hold monitor value
  % handle array/cell of axes
  for index=1:2:numel(varargin) % loop on requested properties
    name = varargin{index}; % axis rank as numeric or string
    value= varargin{index+1};
    if ischar(value) && any(strcmp(value, {'biggest','largest','first','shortest','simplest'}))
      value = findfield(s, '', [ value ' numeric' ]);
    end
    if ~ischar(name) && ~iscellstr(name) && ~isnumeric(name)
      error([ mfilename ': SETAXIS works with axis rank given as char/cellstr/scalar. The ' num2str(index) '-th argument is of type ' class(name) ]);
    end
    if ischar(name), name = cellstr(name);
    elseif isnumeric(name), name = num2cell(name);
    end
    for n_index=1:numel(name)
      % setaxis(s, 'rank'): set the axis value
      % setaxis(a, 'Signal') -> setaxis(a, 0) = Signal/Monitor
      % setaxis(a, 'Error')                   = Error/Monitor
      get_mon = false; sig=[]; err=[];
      % set the alias definition
      if ischar(name{n_index})
        if isscalar(name{n_index})
          if strcmp(name{n_index},'0')  % Signal definition
            s = builtin('subsasgn',s, struct('type','.','subs','Signal'), value);
            s.Private.cache.size = [];
          elseif isfinite(str2num(name{n_index}))
            s.Axes{str2num(name{n_index})}=value;  % set definition in Axes (cell)
          else
            error([ mfilename ': invalid axis rank ''' name{n_index} '''. Should be ''0'' to ''9'' or ''Signal'' or ''Error''.' ]);
          end
        elseif strcmp(name{n_index}, 'Signal')
          name{n_index} = 0;  % then will use rank as number=0, not char
          get_mon = true;
        elseif strcmp(name{n_index}, 'Error')
          get_mon = true; % see below for actual set
        end
      end
      % special case when we need the Monitor value
      if get_mon && isempty(m)
        m = subsref_single(s, 'Monitor'); % follow links -> value
        if ~isnumeric(m), m=1; end
      end
      % second test for 'Error/Monitor' (and now we have Monitor - shared with 'Signal' case)
      if ischar(name{n_index}) && strcmp(name{n_index}, 'Error')
        if isnumeric(value) && (isscalar(m) || isequal(size(m),size(value)))
          value = value.*m;
        end
        s = subsasgn_single(s, 'Error', value); % follow links -> value
      end
      % setaxis(s, rank): set the axis value
      if isnumeric(name{n_index}) && isscalar(name{n_index}) && name{n_index} >= 0
        % set the alias value: interpret result using our subsasgn (follow links)
        if name{n_index} == 0 % {0}=Signal/Monitor
          if isnumeric(value) && (isscalar(m) || isequal(size(m),size(sig)))
            value=value.*m;
          end
          s = subsasgn_single(s, 'Signal', value);
          s.Private.cache.size = [];
        else % Axis
          % set the axis value
          if iscell(s.Axes) && numel(s.Axes) >= name{n_index}
            if ischar(s.Axes{name{n_index}}) && ~ischar(value)
              s = subsasgn_single(s, s.Axes{name{n_index}}, value); % follow link for assignment
            else
              if numel(s.Axes{name{n_index}}) ~= numel(value) % value has changed significantly, request check
                s.Private.cache.check_requested = true;
              end
              s.Axes{name{n_index}} = value; % set definition from char link/alias
            end
          end
        end
      end
    end % for n_index

  end % for index
  history(s, mfilename, varargin{:});
