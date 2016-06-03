function  T = Sab_getT(s)
% Sab_getT: search for a Temperature value in a data set
%
% input:
%   s: any iData object, including S(a,b) ones.
  
  T = [];
  if isfield(s,'temperature'), T=getfield(s,'temperature'); end
  
  if isempty(T)
    f = findfield(s,'Temperature','exact numeric');
    if ~isempty(f), T = s.(f{1}); end
  end
  
  if isempty(T)
    f = findfield(s,'T','exact numeric');
    if ~isempty(f), T = s.(f{1}); end
  end

  if isempty(T)
    f = findfield(s,'Temperature','numeric');
    if ~isempty(f), T = s.(f{1}); end
  end
  
  if isempty(T)
    f = findfield(s,'T','numeric');
    if ~isempty(f), T = s.(f{1}); end
  end
  
  if isvector(T), T = mean(T(:)); else T=[]; end
