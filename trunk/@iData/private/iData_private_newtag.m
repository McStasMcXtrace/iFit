function b = iData_private_newtag(a)
% Data_private_newtag: assigns a new tag and creation date to iData object

b = a;
b.Date = datestr(now);  % new object
tag = tempname;
[dummy, tag] = fileparts(tag);
if length(tag) > 9, tag=tag(1:9); end
b.Tag      = tag;
      
if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),b);
end
