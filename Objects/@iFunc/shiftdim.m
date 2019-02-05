function a=shiftdim(a, dim, parname)
% b=shiftdim(a, rank, shift): shifts the axis of specified rank by a value
%
%   @iFunc/shiftdim function to shift iFunc object axes
%
%   shiftdim(model)
%   shiftdim(model, dim)
%   shiftdim(model, dim, 'offset')
%     The new model has an additional parameter 'offset' which controls the offset 
%     along axis of given rank, e.g. x=x+offset when dim=1. The default rank 
%     (when not given or empty) is 1. The default offset parameter name is e.g. 
%     'Offset_x' (when not given or empty).
%   shiftdim(model, dim, '-name')
%     When the offset parameter name starts with '-' (can be a single '-'), 
%     the parameter is subtracted to the axis, so that this latter is centred
%     e.g. x=x-offset. The default parameter name is then 'Centre_x'.
%
% input:  a:     model (iFunc)
%         rank:  axis rank (scalar). Default is 1 (when empty or not given).
%         shift: name of the offset parameter (string). Default is 1 (when empty 
%                or not given). If the 'shift' name starts with '-', the parameter
%                is used to centre the axis.
% output: b: object (iFunc)
%
% Version: $Date$
% See also  iFunc

if nargin < 2, dim=[]; end
if nargin < 3, parname=''; end

if isempty(dim), dim=1; end
if ~ischar(parname) && ~isempty(parname), return; end

% only when shifting an axis which exists in object
if isnumeric(dim) && dim > a.Dimension && a.Dimension > 0, return; end

centre = false;
if ~isempty(parname) && parname(1) == '-'
  centre = true;
  parname = parname(2:end);
end
ax = 'xyztuvw';
if isnumeric(dim),  ax = ax(dim);
elseif ischar(dim); ax = dim;
else return; end

if isempty(parname) && 
  if centre
    parname = ['Centre_' ax ];
  else
    parname = ['Offset_' ax ];
  end
end

% now we add the code which shifts the axes
a.Parameters{end+1} = parname;
index=numel(a.Parameters);
if centre
  code = [ ax '=' ax '-p(' num2str(index) ')' ];
else
  code = [ ax '=' ax '+p(' num2str(index) ')' ];
end
a = iFunc(plus(code, a)); % add code and check for redundent parameter names
  
if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),a);
end
