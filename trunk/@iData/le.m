function c = le(a,b)
% c = le(a,b) : lower-than-or-equal comparison of iData objects
%
%   @iData/le (<=) comparison operator
%     when comparing two iData objects, the monitor weighting is applied.
%
% input:  a: object or array (iData or numeric)
%         b: object or array (iData or numeric)
% output: c: object or array which Signal is the comparison result (iData)
% ex:     c= a<=1; d=find(a<=b);
%
% Version: $Revision: 1.2 $
% See also iData, iData/find, iData/gt, iData/lt, iData/ge, iData/le, iData/ne, iData/eq

if nargin ==1
	b=[];
end
c = iData_private_binary(a, b, 'le');

