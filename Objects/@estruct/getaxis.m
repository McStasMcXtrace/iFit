function v = getaxis(s,varargin)
% a = getaxis(a,index) get estruct axis value or alias.
%
%   @estruct/getaxis: get estruct axis value or alias.
%   
%     getaxis(a, rank)
%       Get axis value (and follow links).
%
%     getaxis(a, 'rank')
%       Get axis alias.

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
  
  % handle array/cell of axes
  for index=1:numel(varargin) % loop on requested properties
    name = varargin{index}; % axis rank as numeric or string
    if ~ischar(name) && ~iscellstr(name) && ~isscalar(name)
      error([ mfilename ': GETAXIS works with axis rank given as char/cellstr/scalar. The ' num2str(index) '-th argument is of type ' class(name) ]);
    end
    if ischar(name), name = cellstr(name);
    elseif isnumeric(name), name = num2cell(name);
    end
    for n_index=1:numel(name)
      % get the axis definition
      if ischar(name{n_index}) && isscalar(name{n_index})
        v{end+1} = s.Axes{str2num(name{n_index})};  % get definition in Axes (cell)
      elseif isnumeric(name{n_index}) && isscalar(name{n_index})
        if name{n_index} == 0 % Signal/Monitor
          if ~isempty(s.Monitor) && ...
             (isscalar(s.Monitor) || all(size(s.Monitor) == size(s.Signal)))
            v{end+1} = s.Signal./s.Monitor;
          else
            v{end+1} = subsref_single(s, struct('type','.','subs','Signal'));
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
