function s = camproj(a,dim)
% s = camproj(a,dim) : computes the projection of iData objects elements
%
%   @iData/camproj function to compute the projection of the elements of the data set
%     camproj(a,dim) projects along axis of rank dim. All other axes are removed.
%       If dim=0, projection is done on all axes and the total is returned as a scalar value. 
%       camproj(a,1) projects on first dimension (columns).
%       camproj is the complementary to sum.
%
% input:  a: object or array (iData/array of)
%         dim: dimension to project (int)
% output: s: projection of elements (iData/scalar)
% ex:     c=camproj(a);
%
% Version: $Revision: 1.4 $
% See also iData, iData/plus, iData/prod, iData/cumsum, iData/mean, iData/sum

if ~isa(a, 'iData')
  iData_private_error(mfilename,['syntax is sum(iData, dim)']);
end

if nargin < 2, dim=1; end
% handle input iData arrays
if length(a(:)) > 1
  s = {};
  for index=1:length(a(:))
    s{index} = camproj(a(index), dim);
  end
  s = reshape(s, size(a));
  return
end

s=get(a,'Signal');
[link, label]          = getalias(a, 'Signal');
cmd = a.Command;
b=copyobj(a);
setaxis(b, [], getaxis(b)); % delete all axes
if dim == 0
  for index=1:ndims(a)
    s = sum(s, index);
  end
  return  % scalar
else
  % accumulates on all axes except the rank specified
  for index=1:ndims(a)
    if index~=dim, s = sum(s,index); end
  end
  setaxis(b, 1, getaxis(a, num2str(dim)));
  setalias(b,'Signal', s, [ 'projection of ' label ]);     % Store Signal
end

if dim == 1, 
	s=size(b);
	if s(1) == 1, b=transpose(b); end
end

b.Command=cmd;
b = iData_private_history(b, mfilename, b, dim);
s = b;

