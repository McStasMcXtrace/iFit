function a = load_nmoldyn(filename)
% function a=load_nmoldyn(filename)
%
% Returns an iData style dataset from an nMoldyn file
%
% Version: $Date$
% See also: iData/load, iLoad, save, iData/saveas

if ~isa(filename,'iData')
  a = iData(iLoad(filename)); % no post-processing
else
  a = filename; 
end
clear filename

% handle input iData arrays
if numel(a) > 1
  for index=1:numel(a)
    a(index) = feval(mfilename, a(index));
  end
  return
end

% ======================================================================
% searches for some 'known' nMoldyn symbols after import
if isempty(findstr(a, 'nmoldyn')), 
  return; 
end

% 'jobinfo'
f = {'GlobalAttributes.Value','jobinfo','header'};
s = findfield(a, f);
for index=f
  if isfield(a, index{1}), s{end+1} = get(a,index{1}); end
end

% add any other 'fields' as aliases
if ~isempty(s), 
  a=setalias(a, 'jobinfo',s{1},'nMoldyn configuration');
else 
  return;  % no nMoldyn JobInfo stuff
end

% nMoldyn results: last search defines the Signal and Axes: we prefer S(q,w)
a = load_nmoldyn_signal_axes(a, 'pdf_total', 'r');
a = load_nmoldyn_signal_axes(a, 'msd_total', {'times','time'});
a = load_nmoldyn_signal_axes(a, 'dos_total', 'frequency');
a = load_nmoldyn_signal_axes(a, 'atomic_density', {'times','time'});
a = load_nmoldyn_signal_axes(a, 'temperature', {'times','time'});
a = load_nmoldyn_signal_axes(a, 'eisf_total', {'q','k'});
a = load_nmoldyn_signal_axes(a, {'Sq_total','ssf_total'}, {'q','k'});
a = load_nmoldyn_signal_axes(a, {'Fqt_total','f_q_t_total'},{'q','k'}, {'times','time'});
a = load_nmoldyn_signal_axes(a, {'Sqw_total','s_q_f_total'}, {'q','k'}, 'frequency');

% ==============================================================================

function a = load_nmoldyn_signal_axes(a, signal, ax1, ax2)
  % locates items in signal, attach the first match to Signal
  % then locates 'ax1' and 'ax2' and do the same

  if nargin < 3, ax1=''; end
  if nargin < 4, ax2=''; end
  signal = cellstr(signal); ax1= cellstr(ax1); ax2=cellstr(ax2);
  
  % assign Signal looking for member names in object
  [a,s] = load_nmoldyn_search_token(a, signal);
  if isempty(s), return; end
  a = setaxis(a, 0, s);  % assign 0-th axis=Signal to the alias found above
  label(a, 0, s);

  % assign 1st axis looking for member names in object
  [a,s] = load_nmoldyn_search_token(a, ax1);
  if isempty(s), return; end
  a = setaxis(a, 1, s);

  a = setalias(a,'Error',0);

  % assign 2nd axis looking for member names in object
  [a,s] = load_nmoldyn_search_token(a, ax2);
  if isempty(s), return; end
  a = setaxis(a, 2, s);
  
  % check if axes are to be swaped, looking at sizes
  x = getaxis(a,1);
  y = getaxis(a,2);
  if isvector(x) && numel(x) == size(a,2) && isvector(y) && numel(y) == size(a,1)
    x = getaxis(a, '1');
    y = getaxis(a, '2');
    a = setaxis(a, 2, x);
    a = setaxis(a, 1, y);
  end

function [a,alias] = load_nmoldyn_search_token(a, token)
  % searches for a token, and if found checks that a corresponding alias
  % exists, or creates it.
  alias = [];
  if isempty(token), return; end
  
  alias = isfield(a, token);        % use predefined aliases ?
  if any(alias)
    alias = token{find(alias,1)};   % first alias that matches 'token' items
  else
    alias = findfield(a, token, 'exact numeric'); % search in all fields ?
    if ~isempty(alias)
      alias=alias{1};               % first full link e.g. 'Data.<blah>' that matches signal
      a = setalias(a, token{1}, alias);  % a new alias
    end
  end
