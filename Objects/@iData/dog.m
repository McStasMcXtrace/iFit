function s = dog(dim,a, varargin)
% DOG Split object into slices.
%   S=DOG(DIM, A) split objects A into slices along dimension DIM. DIM must be 
%   within dimensionality. Resulting split objects are returned within an array. 
%   This method is the counterpart to CAT.
%
%   S=DOG(DIM, A,B,...) and  S=DOG(DIM, [A,B,...]) split all objects in array.
%
% Example: a=iData(peaks); c=dog(1,a); numel(c) == size(a,1)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/plus, iData/prod, iData/cumcat, iData/mean

if length(varargin) >= 1  % syntax: dog(dim,a,b,c,...)
  s=a;
  for index=1:length(varargin)
    s = [ s ; varargin{index} ];
  end
  clear varargin
  s = dog(dim, s);
  return
end

% syntax is now: dog(dim,[a(:)])
if nargin == 1 & isa(dim, 'iData') & numel(dim) >= 1 % syntax: dog([a])
  s = dog(1, dim);
  return
end

if ~isa(a, 'iData')
  error([mfilename,': syntax is dog(dim, iData, ...)']);
end

% removes warnings during interp


if dim <= 0, dim=1; end
% syntax is now: dog(dim,[a b c ... ])
a=a(:); s = [];
if numel(a) > 1
  s = zeros(iData, numel(a), 1);
  for index=1:numel(a)
    s(index) = dog(dim, a(index));
  end
  return
end

if dim > ndims(a)
  warning([ mfilename, ': Can not extract dim=' num2str(dim) ' slices from object ' this.Tag ' with ndims(a)=' num2str(a) ]);
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
  this_s = history(this_s, mfilename, dim, a);
  s = [ s ; this_s ];
end

