function result=test_iFunc_get

a=gauss;
b=get(a, 'Amplitude');
[signal,a]=feval(a,'guess'); p=a.ParameterValues;
c=get(a, 'Amplitude');
e=get(a,'Expression');
if isempty(b) && isscalar(c) && c == a.Amplitude
  result = [ 'OK     ' mfilename ];
else
  result = [ 'FAILED ' mfilename ];
end
