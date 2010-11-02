function s = dog(dim,a, varargin)
% s= dog(dim,a) : split iData objects elements along dimension
%
%   @iData/dog function to split iData objects elements along dimension dim
%     dog(dim,a) split along axis of rank dim, which must be within dimensionality. 
%     Resulting split objects are returned within an array. This method is the counterpart to cat.
%
% input:  a: object or array (iData)
%         dim: dimension to split (int)
% output: s: split data set (iData array)
% ex:     c=dog(1,a,b); c=dog(1,[ a b ]); 
%
% Version: $Revision: 1.1 $
% See also iData, iData/plus, iData/prod, iData/cumcat, iData/mean

if length(varargin) >= 1  % syntax: dog(dim,a,b,c,...)
  s=a(:);
  for index=1:length(varargin)
    s = [ s ; varargin{index} ];
  end
  s = dog(dim, s);
  return
end

% syntax is now: dog(dim,[a(:)])
if nargin == 1 & isa(dim, 'iData') & length(dim) >= 1 % syntax: dog([a])
  s = dog(1, dim);
  return
end

if ~isa(a, 'iData')
  iData_private_error(mfilename,['syntax is dog(dim, iData, ...)']);
end

% removes warnings during interp
try
  warn.set = warning('off','iData:setaxis');
  warn.get = warning('off','iData:getaxis');
catch
  warn = warning('off');
end

% syntax is now: dog(dim,[a b c ... ])
a=a(:); s = [];
if length(a) > 1
  
  for index=1:length(a)
    s = [ s ; dog(dim, a(index)) ];
  end
  return
end

if dim > ndims(a)
  iData_private_warning(mfilename,[ 'Can not extract dim=' num2str(dim) ' slices from object ' this.Tag ' with ndims(a)=' num2str(a) ]);
  return
end

% prepare the cell of indices to be sent to sub2ind
sub=cell(1,ndims(a));
for index=1:ndims(a)
  if dim ~= index, sub{index}=1:size(a,index); end
end

x = getaxis(a, dim);

for index=1:size(a, dim)
  sub{dim}=index;
  sr.type='()';
  sr.subs=sub;
  this_s=subsref(a,sr);
  setaxis(this_s, dim, [ 'Axis_' num2str(dim) ], x(index));
  this_s = iData_private_history(this_s, mfilename, dim, a);
  s = [ s ; this_s ];
end



% reset warnings during interp
try
  warning(warn.set);
  warning(warn.get);
catch
  warning(warn);
end


