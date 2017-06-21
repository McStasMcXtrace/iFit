function result=test_iFunc_save

a=iFunc('x*p(1)+p(2)'); 
h = save(a);
b = load(iFunc,h);
if ischar(h) && isa(b, 'iFunc') && strcmp(a.Expression, b.Expression)
  result = [ 'OK     ' mfilename ];
  delete(h)
else
  result = [ 'FAILED ' mfilename ];
end
