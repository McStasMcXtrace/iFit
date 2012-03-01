function b = cumsum(a,dim)
% s = cumsum(a,dim) : computes the cumulative sum of iData objects elements
%
%   @iData/cumsum function to compute the cumulative sum of the elements of the data set
%     cumsum(a,dim) accumulates along axis of rank dim.
%
% input:  a: object or array (iData)
%         dim: dimension to accumulate (int)
% output: s: accumulated sum of elements (iData)
% ex:     c=cumsum(a);
%
% Version: $Revision: 1.8 $
% See also iData, iData/plus, iData/sum, iData/prod, iData/cumprod

% handle input iData arrays
if isa(a, 'iData') && numel(a) > 1
  b = [];
  for index=1:numel(a)
    if nargin == 1
      b = [ b feval(mfilename, a(index)) ];
    else
      b = [ b feval(mfilename, a(index), dim) ];
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
  iData_private_error(mfilename,['syntax is cumsum(iData, dim)']);
end
s = iData_private_cleannaninf(get(a,'Signal'));
e = iData_private_cleannaninf(get(a,'Error'));
m = iData_private_cleannaninf(get(a,'Monitor'));

b = setalias(b, 'Signal',   cumsum(s,dim), [ 'cumsum(' sl ','  num2str(dim) ')' ]);
b = setalias(b, 'Error',    cumsum(s+e/2, dim)-cumsum(s-e/2, dim), [ 'cumsum(Error,'  num2str(dim) ')' ]);
b = setalias(b, 'Monitor',  cumsum(m,dim), [ 'cumsum(Monitor,' num2str(dim) ')' ]);
b.Command=cmd;
b = iData_private_history(b, mfilename, a, dim);  

