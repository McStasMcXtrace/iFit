function s = axescheck(s)
% AXESCHECK check the Signal and Axes of an object
%
%   axescheck(s)
%     Check Signal and Axes. When not set, the largest numeric field is assigned
%     as the Signal, and missing Axes are searched accordingly.
%     Monitor and Error are also checked to match Signal dimension.

s.Private.cache.check_requested = 0; % ok, we are working on this

% get the size of the Signal, Monitor, Error (follow links)
signal_sz = size(subsref_single(s, struct('type','.','subs','Signal'))); 
monitor_sz= size(subsref_single(s, struct('type','.','subs','Monitor')));
error_sz  = size(subsref_single(s, struct('type','.','subs','Error')));
axes_id   = [];

% get all numeric fields (e.g. from findfield cache), sorted by decreasing size
[fields, ~, dims, sz] = findfield(s,'','numeric');

% define Signal,Error,Monitor when not yet so
if all(signal_sz == 0) || all(error_sz == 0) || all(monitor_sz == 0) % empty ?
  % identify signal, error, monitor
  [signal_id, error_id, monitor_id, axes_id] = axescheck_find_signal(s, fields, dims, sz);
  % change those which are not set (set definition)
  if all(signal_sz == 0) && ~isempty(signal_id)
    s = builtin('subsasgn', s, struct('type','.','subs','Signal'),  fields{signal_id}); 
    signal_sz  = sz{signal_id};
  end
  if all(error_sz == 0) && ~isempty(error_id)
    s = builtin('subsasgn', s, struct('type','.','subs','Error'),   fields{error_id});
    error_sz   = sz{error_id};
  end
  if all(monitor_id == 0) && ~isempty(monitor_id)
    s = builtin('subsasgn', s, struct('type','.','subs','Monitor'), fields{monitor_id}); 
    monitor_sz = sz{monitor_id};
  end
end

if isempty(signal_sz), return; end

% check Error and Monitor: must be size(Signal) or expression or scalar
if all(error_sz) && prod(signal_sz) ~= prod(error_sz) && ~all(error_sz ==1)
  warning([ mfilename ': invalid Error (size ' mat2str(error_sz) ...
    ' does not match Signal ' mat2str(signal_sz) ') in object ' s.Tag ])
  s = builtin('subsasgn', s, struct('type','.','subs','Error'), []);
end
if all(monitor_sz) && prod(signal_sz) ~= prod(monitor_sz) && ~all(monitor_sz ==1)
  warning([ mfilename ': invalid Monitor (size ' mat2str(monitor_sz) ...
    ' does not match Signal ' mat2str(signal_sz) ') in object ' s.Tag ])
  s = builtin('subsasgn', s, struct('type','.','subs','Monitor'), []);
end

% check Axes: must be size(Signal) or expression or size(rank)
% first search amongst axes_id, then blind search.



% ------------------------------------------------------------------------------
function [signal_id, error_id, monitor_id, axes_id] = axescheck_find_signal(self, fields, dims, sz)
  
  signal_id = []; error_id = []; monitor_id=[]; axes_id = [];

  if      isempty(dims), return; 
  elseif isscalar(dims), signal_id=1; return; 
  end
  fields0 =fields; 
  % we identify named numeric fields, and remove them from our definitions
  for index=1:numel(fields) % decreasing size
    if isempty(fields{index}), continue; end  % skip empty fields/already tagged
    % remove fields that can not be used/cross-referenced
    if any(strcmp(fields{index}, {'Date','ModificationDate','Error','Signal','Monitor'}))
      fields0{index} = ''; dims(index) = 0; sz{index} = 0; fields{index} = '';
      continue; 
    end
    if isempty(error_id) && ~strcmp(fields{index}, 'Error') && ~isempty(strfind(lower(fields{index}), 'err'))
      error_id = index;   % if found, dimension must match signal, or scalar
      fields{index} = '';
    elseif isempty(monitor_id) && ~strcmp(fields{index}, 'Monitor') && ~isempty(strfind(lower(fields{index}), 'mon')) ...
      monitor_id = index; % if found, dimension must match signal, or scalar
      fields{index} = '';
    elseif ~isempty(strfind(lower(fields{index}), 'axis')) ...
      || ~isempty(strfind(lower(fields{index}), 'axes'))
      axes_id(end+1) = index; % if found, dimension must match signal, size(rank)
      fields{index} = '';
    elseif isempty(signal_id) && ~strcmp(fields{index}, 'Signal') ...
      && (~isempty(strfind(lower(fields{index}), 'signa')) ...
      || ~isempty(strfind(lower(fields{index}), 'int')) || ~isempty(strfind(lower(fields{index}), 'amp')))
      signal_id = index;  % signal is named
      fields{index} = '';
    end
  end
  fields = fields0;
  if ~isempty(signal_id) && any(signal_id == [ error_id monitor_id axes_id ])
    signal_id = []; % wrong guess for Signal
  end
  % signal not defined ? get the first biggest field available
  if isempty(signal_id)
    maxdim = dims;
    maxdim([ error_id monitor_id axes_id ]) = []; % remove found named fields
    signal_id = find(maxdim == max(dims),1);  % the first biggest
  end
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

