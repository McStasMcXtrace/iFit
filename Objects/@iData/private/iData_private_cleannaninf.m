function s=iData_private_cleannaninf(s, ratio)
% iData_private_cleannaninf: clean NaNs and Infs from a numerical field
%

if nargin < 2, ratio=10; end
if isa(s,'iData')
  if numel(s) > 1
    a = [];
    for index=1:numel(s)
      a = [ a ; feval(mfilename, s(index), ratio) ];
    end
    s = a;
  else
    s=set(s,'Signal',iData_private_cleannaninf(get(s,'Signal'), ratio));
    s=set(s,'Error', iData_private_cleannaninf(get(s,'Error'), ratio));
  end
  return
end
  
if isnumeric(s)
  S = s(:);
  if all(isfinite(S)), return; end
  index_ok     = find(isfinite(S));

  maxs = max(S(index_ok));
  mins = min(S(index_ok));

  S(isnan(S)) = 0;
  if ~isempty(mins)
    if mins<0, S(find(S == -Inf)) = mins*ratio;
    else       S(find(S == -Inf)) = mins/ratio; end
  end
  if ~isempty(maxs)
    if maxs>0, S(find(S == +Inf)) = maxs*ratio;
    else       S(find(S == +Inf)) = maxs/ratio; end
  end
  s = double(reshape(S, size(s)));
end
