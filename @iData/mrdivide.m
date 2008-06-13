function c = mrdivide(a,b)
% c = mrdivide(a,b) : computes the matrix division of iData objects
%
%   @iData/mrdivide function to compute the matrix division of data sets=a/b
%
% input:  a: object or array (iData or numeric)
%         b: object or array (iData or numeric)
% output: c: object or array (iData)
% ex:     c=a/2;
%
% Version: $Revision: 1.2 $
% See also iData, iData/minus, iData/plus, iData/times, iData/rdivide, iData/power

if isscalar(a) | isscalar(b)
  c = iData_private_binary(a, b, 'rdivide');
else
  iData_private_error(mfilename,[ 'Matrix iData matrix division not supported yet' ]);
end

