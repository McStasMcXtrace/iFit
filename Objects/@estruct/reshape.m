function a = reshape(a, varargin)
% c = reshape(a) : reshape the object Signal
%
%   @estruct/reshape function to reshape the object Signal array
%     reshape(a, m,n,p,...) 
%       reshapes the Signal as an m*n*p*... array the number of elements in
%       the initial Signal must be m*n*p*...
%     reshape(a, [m n p ...]) 
%       is the same thing as above
%
%     The resulting object has the elements of the initial data reordered
%     so that the final size is that requested. The number of elements must
%     not change. To change the number of elements, use estruct/resize or
%     estruct/reducevolume instead. To permute dimensions, use estruct/permute.
%
% input:  a:   object or array (estruct)
%         m,n,p...: dimensions (integers)
% output: c: object or array (estruct)
% ex:     a=estruct(peaks(60)); b=reshape(a, 75, 48);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/squeeze, estruct/size, estruct/permute, estruct/resize,
% estruct/reducevolume

% first get dimensions from varargin
dims = [];
for index=1:length(varargin)
  dims = [ dims varargin{index} ];
end
if isempty(dims), return; end

% handle estruct array: use built-in reshape
if numel(a) > 1
  a = builtin(mfilename, a, dims);
  return
end

% use reshape on Signal, Error, Monitor
sz = size(a);
if prod(sz) ~= prod(dims)
  estruct_private_error(mfilename,[ 'To RESHAPE the number of elements must not change. Object ' ...
      a.Tag ' "' a.Title ' has dimension ' mat2str(size(a)) ' but requested to reshape into ' ...
      mat2str(dims) '. You can rather try the estruct/resize method.' ]) ;
end

a  = unary(a, 'reshape', dims(:)');

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
      x = estruct_private_resize(x, new_sa);
      a = setaxis(a, index, x);
    end
  end
end

if nargout == 0 && ~isempty(inputname(1))
  assignin('caller',inputname(1),a);
end
