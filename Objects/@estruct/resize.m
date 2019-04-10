function b = resize(a, varargin)
% c = resize(a) : resize/rebin the object Signal
%
%   @estruct/resize function to resize/rebin the object Signal array
%     resize(a, m,n,p,...) 
%       resizes the Signal as an m*n*p*... array using a N-D discrete
%       cosine transform.
%     resize(a, [m n p ...])
%       is the same thing as above.
%     resize(a) 
%       interpolates the data using a discrete cosine transform, retaining
%       object dimensions. To further 'clean' the object, you may use
%       estruct/fill, which is slower.
%
% The new dimensions can change the total number of elements, and the
% dimensionality. To reduce uniformly the dimension of the object, you may
% also use estruct/reducevolume. To reshape a data set keeping its
% dimensions, you may use estruct/reshape. To permute dimensions, use
% estruct/permute. To remove singleton dimensions, use estruct/squeeze.
%
% input:  a:   object or array (estruct)
%         m,n,p...: dimensions (integers)
% output: c: object or array (estruct)
% ex:     a=estruct(peaks); b=resize(a, 10,90);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/squeeze, estruct/size, estruct/permute, estruct/reshape,
% estruct/reducevolume, estruct/fill, estruct/resize

% first get dimensions from varargin
dims = [];
for index=1:length(varargin)
  dims = [ dims varargin{index} ];
end
if isempty(dims), dims=size(a); end

% handle estruct array: use built-in reshape
if numel(a) > 1
  b = builtin(mfilename, a, dims);
  return
end

% use reshape on Signal, Error, Monitor
sz = size(a);

% special case of 1D objects resized with higher dimensionality
if ndims(a) == 1 && length(dims) > 1
  b=copyobj(a);
  dims(1:2) = ceil(dims(1:2) ./ size(a)); % new dimension for replication
  set(b, 'Signal', repmat(subsref(a,struct('type','.','subs','Signal')), dims));
else
  b  = unary(copyobj(a), 'estruct_private_resize', dims);
end

% then update axes
myisvector=@(c)length(c) == numel(c);
for index=1:length(dims)
  x  = getaxis(a, index);
  sz = size(x);
  if length(sz) >= index && ndims(a) > 1
    if myisvector(x), x=x(:); sa = [length(x) 1]; else sa = size(x); end
    if ~myisvector(x), new_sa = dims; 
    else             new_sa = [ dims(index) 1]; end

    % resize axis if changed
    if ~isequal(sa, new_sa)
      x = estruct_private_resize(x, new_sa);
      b = setaxis(b, index, x);
    end
  end
end

% and make sure we squeeze object singleton dimensions
b = squeeze(b);

if nargout == 0 && ~isempty(inputname(1))
  assignin('caller',inputname(1),b);
end
