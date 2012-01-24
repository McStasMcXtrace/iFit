function b = cumprod(a,dim)
% s = cumprod(a,dim) : computes the cumulative product of iData objects elements
%
%   @iData/cumprod function to compute the cumulative product of the elements of the data set
%     cumprod(a,dim) operates along axis of rank dim.
%
% input:  a: object or array (iData)
%         dim: dimension to accumulate (int)
% output: s: accumulated product of elements (iData)
% ex:     c=cumprod(a);
%
% Version: $Revision: 1.7 $
% See also iData, iData/plus, iData/sum, iData/prod, iData/cumprod

% handle input iData arrays
if isa(a, 'iData') & length(a(:)) > 1
  b = a(:);
  for index=1:length(a(:))
    if nargin == 1
      b(index) = cumprod(a(index));
    else
      b(index) = cumprod(a(index), dim);
    end
  end
  b = reshape(b, size(a));
  return
end
cmd = a.Command;
b = copyobj(a);
[sn, sl] = getaxis(a, 'Signal');
if nargin == 1
  dim=1;
elseif nargin ~= 2
  iData_private_error(mfilename,['syntax is cumprod(iData, dim)']);
end
s = iData_private_cleannaninf(get(a,'Signal'));
e = iData_private_cleannaninf(get(a,'Error'));
m = iData_private_cleannaninf(get(a,'Monitor'));

b = setalias(b, 'Signal',   cumprod(s,dim), [ 'cumprod(' sl ','  num2str(dim) ')' ]);
b = setalias(b, 'Error',    cumprod(s+e/2,dim)-cumprod(s-e/2,dim), [ 'cumprod(Error,'   num2str(dim) ')' ]);
b = setalias(b, 'Monitor',  cumprod(m,dim), [ 'cumprod(Monitor,' num2str(dim) ')' ]);
b.Command=cmd;
b = iData_private_history(b, mfilename, a, dim);  

