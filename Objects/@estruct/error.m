function e = error(a)
% ERROR  Return the Error bar on the Signal in object

try
  e = subsref(a,struct('type','.','subs','Error'));
catch
  e = [];
end
