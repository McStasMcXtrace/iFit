function c = mtimes(a,b)
% c = mtimes(a,b) : computes the matrix product of iData objects
%
%   @iData/mtimes function to compute the matrix product of data sets=a*b
%      the number of columns of a must equal the number of rows of b.

%
% input:  a: object or array (iData or numeric)
%         b: object or array (iData or numeric)
% output: c: object or array (iData)
% ex:     c=a*2;
%
% Version: $Revision: 1.3 $
% See also iData, iData/minus, iData/plus, iData/times, iData/rdivide, iData/power

if isscalar(a) | isscalar(b)
  c = iData_private_binary(a, b, 'times');
elseif ndims(a) == 2 & ndims(b) == 2
  if size(a,2) ~= size(b,1)
    iData_private_error(mfilename,[ 'the number of columns of a (' num2str(size(a,2)) ') must equal the number of rows of b (' num2str(size(b,1)) ').' ]);
  end
  c = copyobj(a);
  c1 = dot(getaxis(a,2) , getaxis(b,1));
  c2 = dot(getaxis(a,1) , getaxis(b,2));
  c0 = getaxis(a,0)  * getaxis(b,0);
  [a1n, a1l] = getaxis(a, '1');
  [a2n, a2l] = getaxis(a, '2');
  [a0n, a0l] = getaxis(a, '0');
  [b1n, b1l] = getaxis(b, '1');
  [b2n, b2l] = getaxis(b, '2');
  [b0n, b0l] = getaxis(b, '0');
  ae = get(a, 'Error'); am = get(a, 'Monitor');
  be = get(b, 'Error'); bm = get(b, 'Monitor');
  setalias(c,[ a2n '_' b1n ],c1,[ a2l ' ' b1l ]);
  setalias(c,[ a1n '_' b2n ],c2,[ a1l ' ' b2l ]);
  setalias(c, 'Signal', c0, [ a0l '*' b0l ]);
  setalias(c, 'Error',   ae * be);
  setalias(c, 'Monitor', am * bm);
  c = iData_private_history(c, mfilename, a,b);
else
  iData_private_error(mfilename,[ 'Matrix iData multiplication not supported for iData object of dimensions ' num2str(ndims(a)) ' and ' num2str(ndims(b)) ]);
end


