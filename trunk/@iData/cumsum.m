function s = cumsum(a,dim)
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
% Version: $Revision: 1.3 $
% See also iData, iData/plus, iData/sum, iData/prod, iData/cumprod

% handle input iData arrays
if isa(a, 'iData') & length(a(:)) > 1
  s = a(:);
  for index=1:length(a(:))
    if nargin == 1
      s(index) = cumsum(a(index));
    else
      s(index) = cumsum(a(index), dim);
    end
  end
  s = reshape(s, size(a));
  return
end

s = copyobj(a);
[sn, sl] = getaxis(a, '0');
if nargin == 1
  s = setalias(s, 'Signal', cumsum(get(s,'Signal')), [ 'cumsum(' sl ')' ]);
  s = iData_private_history(s, mfilename, a);  
else
  if isa(a, 'iData')
    s = setalias(s, 'Signal', cumsum(get(s,'Signal'), dim), [ 'cumsum(' sl ',' num2str(dim) ')' ]);
    s = iData_private_history(s, mfilename, a, dim);  
  else
    iData_private_error(mfilename,['syntax is cumsum(iData, dim)']);
  end
end

