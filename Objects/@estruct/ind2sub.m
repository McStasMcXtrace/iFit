function data = ind2sub(s, varargin)
% IND2SUB Get indexed array element in a single object.
%   IND2SUB(S,I,J, ...) returns S.Signal(I,J,...)
%   It is equivalent to accessing directly the indexed element in arrays,
%   except when the array is of lenght 1.
%   When S is a single object, S(1) would return S itself,
%   whereas ind2sub(S,1) returns the first element of its 'Signal'.
%
% Example: a=estruct(peaks); isnumeric(ind2sub(a,1))
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/disp, estruct/get, estruct/size
 
  if nargin < 2, return; end

  subs = substruct('()', varargin);
  if length(s) > 1, 
    data=subsref(s, subs);
  else
    if isempty(s), data=[]; 
    else 
    	s = get(s,'Signal');
    	data = subsref(s, subs);
    end
  end    
