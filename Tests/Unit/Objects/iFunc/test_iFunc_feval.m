function result=test_iFunc_feval

a=gauss;
[s,a,ax] = feval(a,NaN);
s_guess  = feval(a,'guess');
if length(a.ParameterValues) == length(a.Parameters) ...
&& length(s) == length(ax{1})
  result = [ 'OK     ' mfilename ];
else
  result = [ 'FAILED ' mfilename ];
end
