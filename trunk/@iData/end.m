function b = end(s,k,n)
% b = end(s,index,n) : end value for iData objects
%
%   @iData/end function defines end value for iData
%   returns the length of rank 'index' among total dimensions 'n' in object 's'.
%
% Version: $Revision: 1.3 $
% See also iData

% EF 27/07/00 creation
% EF 23/09/07 iData implementation

if length(s(:)) > 1
  if n == 1, b=length(s(:)); else b=size(s,k); end
  return
end

if length(size(s.Signal)) < n
  iData_private_error(mfilename, ['input iData object ' inputname(1) ' ' b.Tag ' has a size [' num2str(size(s)) '] but the dimension ' n ' is requested.' ]);
end

b = size(s.Signal,k);
  
