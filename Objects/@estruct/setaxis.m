function s = setaxis(s,varargin)
% a = setaxis(a,index) set estruct axis value or alias.
%
%   @estruct/setaxis: set estruct axis value or alias.
%   
%     setaxis(a, rank, value)
%       Set axis value (and follow links).
%
%     setaxis(a, 'rank', value)
%       Set axis alias/link/definition.
%
%     setaxis(a)
%       Check signal and axes.

%   An alias is a string/char which allows to link to internal or 
%   external links, as well as evaluated expression, with the following syntax cases:
%     'field'                           a simple link to an other property 'field'
%     'field1.field2...'                a nested link to an other property
%     'file://some_file_path'           a local file URL
%     'http://some_distant_resource'    an HTTP URL (proxy settings may have to be set)
%     'https://some_distant_resource'   an HTTPS URL (proxy settings may have to be set)
%     'ftp://some_distant_resource'     an FTP URL (proxy settings may have to be set)
%     'matlab: some_expression'         some code to evaluate. 'this' refers to the object itself
%
%   File and URL can refer to compressed resources (zip, gz, tar, Z) which are 
%   extracted on-the-fly. In case the URL/file resource contains 'sections', a 
%   search token can be specified with syntax such as 'file://filename#token'.
%
% Version: $Date$

  if nargin == 1, s = axescheck(s); return; end
  
  % handle array of struct
  v = {};
  if numel(s) > 1
    for index=1:numel(s)
      s(index) = setaxis(s(index), varargin{:});
    end
    v = reshape(v, size(s));
    return
  end
  
  % handle array/cell of axes
  for index=1:2:numel(varargin) % loop on requested properties
    name = varargin{index}; % axis rank as numeric or string
    value= varargin{index+1};
    if ~ischar(name) && ~iscellstr(name) && ~isscalar(name)
      error([ mfilename ': SETAXIS works with axis rank given as char/cellstr/scalar. The ' num2str(index) '-th argument is of type ' class(name) ]);
    end
    if ischar(name), name = cellstr(name);
    elseif isnumeric(name), name = num2cell(name);
    end
    for n_index=1:numel(name)
      
      % set the alias definition
      if ischar(name{n_index}) && (strcmpi(name{n_index}, 'Signal') || isscalar(name{n_index}))
        if str2num(name{n_index}) > 0       % Axes
          s.Axes{str2num(name{n_index})} = value; % set definition in Axes (cell)
        elseif strcmpi(name{n_index}, 'Signal') || str2num(name{n_index}) == 0  % Signal
          s = builtin('subsasgn',s, struct('type','.','subs','Signal'), value);
        end
      elseif isnumeric(name{n_index}) && isscalar(name{n_index})
        % set the alias value: interpret result using our subsasgn (follow links)
        if name{n_index} == 0 % {0}=Signal/Monitor
          if ~isempty(s.Monitor) && ...
             (isscalar(s.Monitor) || all(size(s.Monitor) == size(value)))
            s.Signal = value.*s.Monitor;
          else
            s.Signal = value;
          end
        elseif name{n_index}>0 % Axis
          % set the axis value
          s.Axes{name{n_index}} = value;
        end
      end
    end % for n_index
    
  end % for index
  s.Private.cache.check_requested = true;
