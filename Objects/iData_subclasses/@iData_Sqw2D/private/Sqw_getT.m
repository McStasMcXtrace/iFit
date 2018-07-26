function  T = Sqw_getT(s, prop, raw)
% Sqw_getT: search for a property value in a data set
%
%   T = Sqw_getT(s, prop, raw)
%
% input:
%   s: any iData object, including S(q,w) and DOS ones.
%   prop: a list of equivalent properties (numeric) to search for.
%   raw: optional, when specified, use raw output, else compute the mean value.
%
% output:
%   T: the property value (mean value when 'raw' not specified)

  if nargin < 3, raw = false; else raw = true; end
  
  T = [];
  if nargin < 2,    prop = []; end
  if isempty(prop), prop={'Temperature','T'}; end
  if ~iscellstr(prop), prop = { prop }; end
  
  % handle arrays
  if numel(s) > 1
    for index=1:numel(s)
      t = Sqw_getT(s(index), prop);
      if isempty(t), t=nan; end
      T = [ T t ];
    end
    return
  end
  
  for field=prop
    % does it already exist as a property ?
    if isfield(s,field{1}), T=get(s,field{1}); end
    
    % search case sensitive, the first match only
    if isempty(T) || all(T(:)<0) || ~isnumeric(T) || ~isvector(T)
      f = findfield(s,field{1},'exact numeric first');
      if ~isempty(f), T = get(s,f); end
    end
    
    % not found ? search non case sensitive, all matches
    if isempty(T) || all(T(:)<0) || ~isnumeric(T) || ~isvector(T)
      f = findfield(s,field{1},'numeric cache');
      for index=1:numel(f)
        try
          T = get(s,f{1});
          if ~isempty(T) && isnumeric(T) && isvector(T) && all(T(:)>0), break; end
        end
      end
    end

    if ~isempty(T) && isnumeric(T) && isvector(T) && all(T(:)>0), break; end
  end
  
  if raw
    if isvector(T), T = mean(T(:)); else T=[]; end % also for scalars
  end
  
 
