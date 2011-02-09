function s=iData_private_cleannaninf(s)
% iData_private_cleannaninf: clean NaNs and Infs from a numerical field
%

if isnumeric(s)
  S = s(:);
  if all(isfinite(S)), return; end
  index_ok     = find(isfinite(S));

  maxs = max(S(index_ok));
  mins = min(S(index_ok));

  S(isnan(S)) = 0;
  if mins<0, S(find(S == -Inf)) = mins*100;
  else       S(find(S == -Inf)) = mins/100; end
  if maxs>0, S(find(S == +Inf)) = maxs*100;
  else       S(find(S == +Inf)) = maxs/100; end

  s = reshape(S, size(s));
end
