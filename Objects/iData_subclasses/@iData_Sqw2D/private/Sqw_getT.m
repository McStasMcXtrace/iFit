function  T = Sqw_getT(s, prop)
% Sqw_getT: search for a property value in a data set
%
% input:
%   s: any iData object, including S(q,w) and DOS ones.
  
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
    if isfield(s,field{1}), T=get(s,field{1}); end
    
    if isempty(T) || all(T(:)<0) || ~isnumeric(T) || ~isvector(T)
      f = findfield(s,field{1},'exact numeric');
      if ~isempty(f), T = get(s,f{1}); end
    end

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
  
  if isvector(T), T = mean(T(:)); else T=[]; end
  
 
