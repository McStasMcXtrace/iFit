function s = axescheck(s, option)
% AXESCHECK Check the Signal and Axes of an object.
%   AXESCHECK(s) Check Signal and Axes. When not set, the largest numeric
%   field is assigned as the Signal, and missing Axes are searched accordingly.
%   Monitor and Error are also checked to match Signal dimension.
%
%   Usually, the objects are checked automatically for integrity, and it
%   is not required to manually call this method.
%
%   AXESCHECK(s, 'nocheck') makes no check and prevents object to be checked
%   until a further assignment is done. This can be used to assign many properties
%   without testing the object. Only a final assignment will trigger a check.
%
% Example: s = estruct(1:10); isa(axescheck(s),'estruct')
% Version: $Date$ $Version$ $Author$
% see also estruct, getaxis, setaxis

% handle array of struct
if numel(s) > 1
  for index=1:numel(s)
    s(index) = axescheck(s(index));
  end
  return
end

notify(s, 'ObjectUpdated');
s.Private.cache.check_requested = false; % ok, we are working on this

if nargin == 2 && strcmp(option, 'nocheck')
  return
end

%% Check Signal, Monitor, Error ================================================

% get the size of the Signal, Monitor, Error (follow links)
signal_sz = size(s);
monitor_sz= size(subsref_single(s, 'Monitor'));
error_sz  = size(subsref_single(s, 'Error'));
axes_id   = [];

% get all numeric fields (e.g. from findfield cache), sorted by decreasing size
[fields, ~, dims, sz] = findfield(s,'','numeric');
if s.verbose > 1
  for index=1:numel(fields)
    if dims(index) > 1
      disp([ mfilename ': ' num2str(index) ' as ' fields{index} ' numel [' num2str(dims(index)) ']' ])
    end
  end
  disp([ mfilename ': scalars as ' sprintf('%s ', fields{dims == 1}) ]);
  disp([ mfilename ': empty   as ' sprintf('%s ', fields{dims == 0}) ]);
end

% define Signal,Error,Monitor when not yet so
if all(signal_sz == 0) || all(error_sz == 0) || all(monitor_sz == 0) % empty ?
  % identify signal, error, monitor and axes based on names. Signal can also be set as the biggest field.
  [signal_id, error_id, monitor_id, axes_id] = axescheck_find_signal(s, fields, dims, sz);
  if s.verbose > 1
    disp([ mfilename ': called axescheck_find_signal' ]);
    disp([ '  signal:  ' num2str(signal_id) ' as ' fields{signal_id} ])
    if ~isempty(error_id),   disp([ '  error:   ' num2str(error_id)    ' as ' fields{error_id} ]); end
    if ~isempty(monitor_id), disp([ '  monitor: ' num2str(monitor_id)  ' as ' fields{monitor_id} ]); end
    if ~isempty(axes_id),    disp([ '  axes:    ' num2str(axes_id) ]); end
  end
  % change those which are not set (set definition)
  if all(signal_sz == 0) && ~isempty(signal_id)
    s = builtin('subsasgn', s, struct('type','.','subs','Signal'),  fields{signal_id});
    signal_sz  = sz{signal_id};
    if s.verbose > 1,
      disp([ mfilename ': setting Signal = ' fields{signal_id} ' [' num2str(dims(signal_id)) ']' ]);
    end
  end
  if ~isempty(signal_sz)
    if all(error_sz == 0) && ~isempty(error_id)
      s = builtin('subsasgn', s, struct('type','.','subs','Error'),   fields{error_id});
      error_sz   = sz{error_id};
      if s.verbose > 1,
        disp([ mfilename ': setting Error = ' fields{error_id} ' [' num2str(dims(error_id)) ']' ]);
      end
    end
    if all(monitor_id == 0) && ~isempty(monitor_id)
      s = builtin('subsasgn', s, struct('type','.','subs','Monitor'), fields{monitor_id});
      monitor_sz = sz{monitor_id};
      if s.verbose > 1,
        disp([ mfilename ': setting Monitor = ' fields{monitor_id} ' [' num2str(dims(monitor_id)) ']' ]);
      end
    end
  end
end

if isempty(signal_sz), return; end

% check Error and Monitor: must be size(Signal) or expression or scalar
if all(error_sz) && prod(signal_sz) ~= prod(error_sz) && ~all(error_sz ==1)
  if s.verbose
    warning([ mfilename ': invalid Error (size ' mat2str(error_sz) ...
      ' does not match Signal ' mat2str(signal_sz) ') in object ' s.Tag ])
  end
  s = builtin('subsasgn', s, struct('type','.','subs','Error'), []);
end
if all(monitor_sz) && prod(signal_sz) ~= prod(monitor_sz) && ~all(monitor_sz ==1)
  if s.verbose
    warning([ mfilename ': invalid Monitor (size ' mat2str(monitor_sz) ...
      ' does not match Signal ' mat2str(signal_sz) ') in object ' s.Tag ])
  end
  s = builtin('subsasgn', s, struct('type','.','subs','Monitor'), []);
end
s.Private.cache.size = signal_sz; % for faster size execution

%% check Axes: must be size(Signal) or expression or size(rank)
% first search amongst axes_id, then blind search.
ok = axescheck_find_axes(s, fields(axes_id), dims(axes_id), sz(axes_id));
if ok < numel(signal_sz)
  % some axes were not found from reduced Axes search. Try with full object.
  axescheck_find_axes(s, fields,          dims,          sz);
end

% ------------------------------------------------------------------------------
function [signal_id, error_id, monitor_id, axes_id] = axescheck_find_signal(self, fields, dims, sz)
  % AXESCHECK_FIND_SIGNAL Search for Signal, Error, Monitor and Axes based on their names
  %   Signal can also be set as the biggest field.
  %   Return the index in the field cell.
  signal_id = []; error_id = []; monitor_id=[]; axes_id = [];

  if      isempty(dims), return;
  elseif isscalar(dims), signal_id=1; return;
  end

  % clean fields with 'Protected' stuff.
  for index=1:numel(fields) % decreasing size
    if isempty(fields{index}), continue; end  % skip empty fields/already tagged
    % remove fields that can not be used/cross-referenced
    field_root = strtok(fields{index},'.');
    if any(strcmp(fields{index}, self.properties_Base)) ...
    || any(strcmp(fields{index}, self.properties_Protected)) ...
    || any(strcmp(fields{index}, {'Error','Signal','Monitor'})) ...
    || any(strcmp(field_root,    {'Attributes','Headers','Labels','postprocess','verbose'}))
      fields{index} = ''; dims(index) = 0; sz{index} = 0;
    end
  end

  % we search for Monitor, Error and Axes by name.
  % we identify named numeric fields, and remove them from our definitions
  for index=1:numel(fields) % decreasing size
    if isempty(fields{index}), continue; end  % skip empty fields/already tagged
    if isempty(error_id) && ~isempty(strfind(lower(fields{index}), 'err'))
      error_id = index;   % if found, dimension must match signal, or scalar
      fields{index} = '';
    elseif isempty(monitor_id) && ~isempty(strfind(lower(fields{index}), 'mon')) ...
      monitor_id = index; % if found, dimension must match signal, or scalar
      fields{index} = '';
    elseif ~isempty(strfind(lower(fields{index}), 'axis')) ...
      ||   ~isempty(strfind(lower(fields{index}), 'axes'))
      axes_id(end+1) = index; % if found, dimension must match signal, size(rank)
      fields{index} = '';
    end
  end

  % We search for the Signal. Its size must be max(dims). Its name must not match
  % 'error','monitor' or any other alias.
  maxdim = max(dims); signal_id_candidate = [];
  for index=1:numel(fields) % decreasing size
    if isempty(fields{index}), continue; end  % skip empty fields/already tagged
    if isempty(signal_id) && ~strcmp(fields{index}, 'Signal') ...
      && (~isempty(strfind(lower(fields{index}), 'signa')) ...
        || ~isempty(strfind(lower(fields{index}), 'int')) ...
        || ~isempty(strfind(lower(fields{index}), 'amp')))
      signal_id_candidate(end+1) = index;  % signal is named
      fields{index} = '';
    elseif dims(index) == maxdim
      signal_id_candidate(end+1) = index;
    else break; % not max size, not a nice name.
    end
  end
  if ~isempty(signal_id_candidate), signal_id = signal_id_candidate(1); end

  % no Signal found: then nothing can be defined, return
  if isempty(signal_id)
    error_id = []; monitor_id=[]; axes_id = [];
    return
  end

  % now check dimensions wrt signal
  if ~isempty(error_id) && ~all(dims(signal_id) == dims(error_id)) && dims(error_id) ~= 1
    error_id = [];  % found error does not match signal nor scalar
  end
  if ~isempty(monitor_id) && ~all(dims(signal_id) == dims(monitor_id)) && dims(monitor_id) ~= 1
    monitor_id = [];  % found monitor does not match signal nor scalar
  end

% ------------------------------------------------------------------------------
function ok = axescheck_find_axes(self, fields, dims, sz)

  ok = 0;
  if isempty(dims), return; end

  % scan Axes/dimensions
  for index=1:ndims(self)
    % get current axis definition
    if index <= numel(self.Axes)
      def = self.Axes{index};
    else def=[]; end

    % when not empty, we check that it is:
    %   size(Signal) or size(rank)
    if ~isempty(def)
      if ischar(def)
        try
          ax_sz = size(get(self, def));
        catch
          ax_sz = [];
        end
      elseif isnumeric(def)
        ax_sz = size(def);
      end

      if axescheck_size_axes(self, index, ax_sz)
        if self.verbose > 1
          disp([ mfilename ': axescheck_find_axes: checked axis[' num2str(index) '] ' num2str(def) ' with size ' mat2str(ax_sz) ]);
        end
        ok = ok + 1;
        continue
      end

      % if we reach here, axis is invalid
      def=[]; self.Axes{index} = [];
    end

    if isempty(def)
    % when empty we search for one amongst fields that is:
    %   size(Signal) or size(rank)

      % get definition of Signal, Error and Monitor, so that we do not use these.
      not_use = self.Axes; % will not use already defined axes
      for findex={'Signal','Error','Monitor'}
        not_use{end+1} = getalias(self, findex{1});
        not_use{end+1} = findex{1};
      end
      not_use = not_use(~cellfun(@isempty, not_use)); % clean empty elements

      for findex = 1:numel(fields)
        % check for field
        if ~any(strcmp(fields{findex}, not_use)) && axescheck_size_axes(self, index, sz{findex})
          % found a valid axis definition/value
          if self.verbose > 1
            disp([ mfilename ': axescheck_find_axes: found axis[' num2str(index) '] as ' fields{findex} ' with size ' mat2str(sz{findex}) ]);
          end
          self.Axes{index} = fields{findex};
          ok = ok + 1;
          break
        end
      end
    end
  end % for

  % ----------------------------------------------------------------------------
  function tf = axescheck_size_axes(self, index, ax_sz)
    tf = false;
    % size(Signal) or size(rank) as vector
    %   NOTE: some data sets have 'central' axes values with numel = dim(Signal)-1
    if (numel(ax_sz) == numel(size(self)) && all(ax_sz == size(self))) ...
    || (numel(ax_sz) >= index && isscalar(find(ax_sz>1)) && prod(ax_sz) == size(self, index)) ...
    || (numel(ax_sz) >= index && isscalar(find(ax_sz>1)) && prod(ax_sz) == size(self, index)-1) ...
      tf = true;  % valid axis definition/value
    end
