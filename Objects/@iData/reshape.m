function b = reshape(a, varargin)
% RESHAPE Reshape object (reorder).
%   B = RESHAPE(A, M,N,P,...) returns the object A as an M*N*P*... array. 
%   An error results if X does not have M*N elements.
%   The resulting object has the elements of the initial data reordered
%   so that the final size is that requested. The number of elements must
%   not change. To change the number of elements, use RESIZE (with interpolation)
%   or REDUCEVOLUME instead. To permute dimensions, use PERMUTE.
%
%   B = RESHAPE(A, [m n p ...]) is the same thing as above.
%
% Example: a=iData(peaks(60)); b=reshape(a, 75, 48); all(size(b) == [75 48])
% Version: $Date$ $Version$ $Author$
% See also iData, iData/squeeze, iData/size, iData/permute, iData/resize,
% iData/reducevolume

% first get dimensions from varargin
dims = []; b = [];
for index=1:length(varargin)
  dims = [ dims varargin{index} ];
end
if isempty(dims), return; end

% handle iData array: use built-in reshape
if numel(a) > 1
  b = builtin(mfilename, a, dims);
  return
end

% use reshape on Signal, Error, Monitor
sz = size(a);
if prod(sz) ~= prod(dims)
  error([ mfilename ': To RESHAPE the number of elements must not change. Object ' ...
      a.Tag ' "' a.Name ' has dimension ' mat2str(size(a)) ' but requested to reshape into ' ...
      mat2str(dims) '. You may rather try the RESIZE method.' ]) ;
end

b  = unary(a, 'reshape', dims(:)');

% then update axes
myisvector = @(c)length(c) == numel(c);
for index=1:length(dims)
  if length(sz) >= index && sz(index) ~= dims(index)
    x  = getaxis(a, index);
    if myisvector(x),
        x=x(:); sa = [length(x) 1];
        new_sa = [ dims(index) 1 ];
    else
        sa = size(x);
        new_sa = dims;
    end

    % resize axis if changed
    if ~isequal(sa, new_sa)
      x = private_resize(x, new_sa);
      a = setaxis(a, index, x);
    end
  end
end
