function result=test_iFunc_feval

a=gauss;
[s,a,ax] = feval(a,NaN);
if length(feval(a,'guess')) == length(a.Parameters) ...
&& length(s) == length(ax{1})
  result = [ 'OK     ' mfilename ];
else
  result = [ 'FAILED ' mfilename ];
end
