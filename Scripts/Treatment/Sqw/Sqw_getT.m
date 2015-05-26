function  T = Sqw_getT(s)
  % Sqw_getT: search for a Temperature value is data set
  
  T = [];
  f = findfield(s,'Temperature','exact numeric');
  if ~isempty(f), T = s.(f{1}); return; end
  
  f = findfield(s,'T','exact numeric');
  if ~isempty(f), T = s.(f{1}); return; end

  f = findfield(s,'Temperature','numeric');
  if ~isempty(f), T = s.(f{1}); return; end
  
  f = findfield(s,'T','numeric');
  if ~isempty(f), T = s.(f{1}); return; end
