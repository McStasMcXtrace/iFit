function [v,lab] = getaxis(s,varargin)
% GETAXIS Get axis definition or value in object.
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
% Version: $Date$ $Version$ $Author$
% See also estruct, fieldnames, findfield, isfield, set, get, getalias, setalias,
%   getaxis, setaxis

  if nargin == 1, v = s.Axes; if nargout > 1, lab = label(s, 1:ndims(s)); end; return; end

  % handle array of struct
  v = {}; lab={};
  if numel(s) > 1
    for index=1:numel(s)
      [v{end+1},lab{end+1}] = getaxis(s(index), varargin{:});
    end
    v   = reshape(v,   size(s));
    lab = reshape(lab, size(s));
    return
  end

  % check object when we evaluate/get some data out of it, and changes were marked.
  if isfield(s.Private,'cache') && isfield(s.Private.cache,'check_requested') ...
    && s.Private.cache.check_requested
    axescheck(s);
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
      if s.verbose > 2
        if isnumeric(name{n_index})
          disp([ mfilename ': DEBUG: get axis ' num2str(name{n_index}) ])
        else
          disp([ mfilename ': DEBUG: get axis ' class(name{n_index}) ' ' char(name{n_index}) ])
        end
      end
      if ischar(name{n_index})
        if isscalar(name{n_index})
          if strcmp(name{n_index},'0')  % Signal definition
            v{end+1} = builtin('subsref',s, struct('type','.','subs','Signal'));

          elseif isfinite(str2num(name{n_index})) && str2num(name{n_index}) <= numel(s.Axes)
            v{end+1} = s.Axes{str2num(name{n_index})};  % get definition in Axes (cell)
          else
            v{end+1} = []; % axis rank beyond dimension of object Signal
          end
          if nargout > 1, lab{end+1} = label(s, str2num(name{n_index})); end
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
        if isempty(m) || ~isnumeric(m) || all(~isfinite(m(:))) || all(~m(:)), m=1; end
      end
      % second test for 'Error/Monitor' (and now we have Monitor - shared with 'Signal' case)
      if ischar(name{n_index}) && strcmp(name{n_index}, 'Error')
        err = subsref_single(s, 'Error'); % follow links -> value
        if isscalar(m) || isequal(size(m),size(err)), v{end+1} = err./m;
        else                                          v{end+1} = err;
        end
        continue
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
          else % axis rank beyond dimension of object Signal
            v{end+1} = [];
          end
          if nargout > 1, lab{end+1} = label(s, name{n_index}); end
          if (isempty(v{end}) || (isscalar(v{end}) && isnan(v{end}))) && name{n_index} <= ndims(s)
            if numel(find(size(s)>1)) == 1
              v{end} = 1:prod(size(s));
            else
              v{end} = 1:size(s, name{n_index});
            end
          end
          % check axis orientation
          n = size(v{end});
          if isscalar(find(n > 1))
            if length(find(size(s) > 1)) ~= 1
              z = ones(1, length(n));
              z(name{n_index}) = max(n);
              if prod(size(v{end})) == prod(z), v{end}   = reshape(v{end}, z); end
            else
              if prod(size(v{end})) == prod(size(s)), v{end} = reshape(v{end}, size(s)); end
            end
          end
          if isnumeric(v{end}) && ~isfloat(v{end})
            v{end} = double(v{end});
          end
        end
      end
    end
  end
  if numel(v) == 1, v = v{1}; end
  if nargout > 1
    if numel(lab) == 1, lab=lab{1}; end
  end
