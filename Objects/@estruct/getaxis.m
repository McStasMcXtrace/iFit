function v = getaxis(s,varargin)
% GETAXIS get axis definition or value in object.
%   GETAXIS(a, rank) Get axis value (and follow aliases). This is equivalent to
%   the syntax a{rank}. The axis rank 0 corresponds with the Signal/Monitor value.
%   The axis of rank 1 corresponds with rows, 2 with columns, 3 with pages, etc.
%
%   GETAXIS(a, 'Signal') Get the Signal/Monitor value, equivalent to a{0} and 
%   getaxis(a, 0).
%
%   GETAXIS(a, 'Error')  Get the Error/Monitor value.
%
%   GETAXIS(a, 'rank') Get axis definition (alias). The axis of rank 0
%   corresponds with the Signal definition.
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
% Example: s=estruct('x',1:10,'y',1:20, 'data',rand(10,20)); setaxis(s,1,'x'); isnumeric(getaxis(s,1))
% Version: $Date$ (c) E.Farhi. License: EUPL.
% See also estruct, fieldnames, findfield, isfield, set, get, getalias, setalias, 
%   getaxis, setaxis

  if nargin == 1, v = getaxis(s, 1:numel(s.Axes)); return; end
  
  % handle array of struct
  v = {};
  if numel(s) > 1
    for index=1:numel(s)
      v{end+1} = getaxis(s(index), varargin{:});
    end
    v = reshape(v, size(s));
    return
  end
  
  m = []; % will hold monitor value
  % handle array/cell of axes
  for index=1:numel(varargin) % loop on requested properties
    name = varargin{index}; % axis rank as numeric or string
    if ~ischar(name) && ~iscellstr(name) && ~isnumeric(name)
      error([ mfilename ': GETAXIS works with axis rank given as char/cellstr/scalar. The ' num2str(index) '-th argument is of type ' class(name) ]);
    end
    if ischar(name), name = cellstr(name);
    elseif isnumeric(name), name = num2cell(name);
    end
    for n_index=1:numel(name)
      % getaxis(s, 'rank'): get the axis value
      % getaxis(a, 'Signal') -> getaxis(a, 0) = Signal/Monitor
      % getaxis(a, 'Error')                   = Error/Monitor
      get_mon = false; sig=[]; err=[];
      % get the alias definition
      if ischar(name{n_index})
        if isscalar(name{n_index})
          if strcmp(name{n_index},'0')  % Signal definition
            v{end+1} = builtin('subsref',s, struct('type','.','subs','Signal'));
          elseif isfinite(str2num(name{n_index}))
            v{end+1} = s.Axes{str2num(name{n_index})};  % get definition in Axes (cell)
          else 
            error([ mfilename ': invalid axis rank ''' name{n_index} '''. Should be ''0'' to ''9'' or ''Signal'' or ''Error''.' ]);
          end
        elseif strcmp(name{n_index}, 'Signal')
          name{n_index} = 0;  % then will use rank as number=0, not char
          get_mon = true;
        elseif strcmp(name{n_index}, 'Error')
          get_mon = true; % see below for actual get
        end
      end
      % special case when we need the Monitor value
      if get_mon && isempty(m)
        m = subsref_single(s, 'Monitor'); % follow links -> value
        if ~isnumeric(m), m=1; end
      end
      % second test for 'Error/Monitor' (and now we have Monitor - shared with 'Signal' case)
      if ischar(name{n_index}) && strcmp(name{n_index}, 'Error')
        err = subsref_single(s, 'Error'); % follow links -> value
        if isscalar(m) || isequal(size(m),size(err)), v{end+1} = err./m;
        else                                          v{end+1} = err;
        end
      end
      % getaxis(s, rank): get the axis value
      if isnumeric(name{n_index}) && isscalar(name{n_index}) && name{n_index} >= 0
        if name{n_index} == 0 % Signal/Monitor
          sig = subsref_single(s, 'Signal');
          if isscalar(m) || isequal(size(m),size(sig)), v{end+1} = sig./m;
          else                                          v{end+1} = sig;
          end
        else % Axis
          % get the axis value
          if iscell(s.Axes) && numel(s.Axes) >= name{n_index}
            def = s.Axes{name{n_index}};
            if ischar(def)
              v{end+1} = subsref_single(s, def);
            else
              v{end+1} = def;
            end
          else % invalid Axis (not assigned): use a vector matching Signal dimension
            v{end+1} = 1:size(s.Signal, name{n_index});
          end
        end
      end
    end
  end
  if numel(v) == 1, v = v{1}; end
