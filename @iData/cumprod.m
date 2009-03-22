function s = cumprod(a,dim)
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
% Version: $Revision: 1.4 $
% See also iData, iData/plus, iData/sum, iData/prod, iData/cumsum

% handle input iData arrays
if isa(a, 'iData') & length(a(:)) > 1
  s = a(:);
  for index=1:length(a(:))
    if nargin == 1
      s(index) = cumprod(a(index));
    else
      s(index) = cumprod(a(index), dim);
    end
  end
  s = reshape(s, size(a));
  return
end
cmd =a.Command;
s = copyobj(a);
[sn, sl] = getaxis(a, '0');
if nargin == 1
  s = setalias(s, 'Signal', cumprod(get(s,'Signal')), [ 'cumprod(' sl ')' ]);
  s.Command=cmd;
  s = iData_private_history(s, mfilename, a);  
else
  if isa(a, 'iData')
    s = setalias(s, 'Signal', cumprod(get(s,'Signal'), dim), [ 'cumprod(' sl ',' num2str(dim) ')' ]);
    s.Command=cmd;
    s = iData_private_history(s, mfilename, a, dim);  
  else
    iData_private_error(mfilename,['syntax is cumprod(iData, dim)']);
  end
end

