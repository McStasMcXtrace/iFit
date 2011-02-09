function b = iData_private_newtag(a)
% Data_private_newtag: assigns a new tag and creation date to iData object

persistent id

b = a;
b.Date = datestr(now);  % new object

% create new tag
if ~exist('id'),  id=0; end
if isempty(id),   id=0; end
if id <=0, 
  id = clock;
  id = fix(id(6)*999999); 
else 
  id=id+1;
  if id > 1e7, id=0; end
end

b.Tag      = [ 'id' sprintf('%0.f', id) ];;
      
if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),b);
end
