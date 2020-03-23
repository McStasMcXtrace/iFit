function e = error(a)
% ERROR  Return the Error bar on the Signal in object.
%   ERROR(s) is similar to s.Error
%
% Example: s=estruct([1 1 1 1]); all(error(s) == 1)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/floor, estruct/ceil, estruct/round

try
  e = subsref(a,struct('type','.','subs','Error'));
catch
  e = [];
end
