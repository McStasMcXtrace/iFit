function  T = Sab_getT(s)
% Sab_getT: search for a Temperature value in a data set
%
% input:
%   s: any iData object, including S(a,b) ones.
  
  T = [];
  
  % handle arrays
  if numel(s) > 1
    for index=1:numel(s)
      t = Sqw_getT(s(index));
      if isempty(t), t=nan; end
      T = [ T t ];
    end
    return
  end
  
  if isfield(s,'temperature'), T=getfield(s,'temperature'); end
  
  if isempty(T)
    f = findfield(s,{'Temperature','T'},'exact numeric');
    if ~isempty(f), T = s.(f{1}); end
  end

  if isempty(T)
    f = findfield(s,{'Temperature','T'},'numeric cache');
    if ~isempty(f), T = s.(f{1}); end
  end
  
  if isvector(T), T = mean(T(:)); else T=[]; end
