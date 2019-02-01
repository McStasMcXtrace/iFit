function a=shiftdim(a, varargin)
% b=shiftdim(a, rank, shift): shifts the axis of specified rank by a value
%
%   @iData/shiftdim function to shift iData object axes
%     This is equivalent to a{rank} = a{rank}+shift;
%
% input:  a:     object or array (iData)
%         rank:  axis rank (scalar)
%         shift: value to shift the axis with (scalar)
% output: b: object or array (iData)
%
% Version: $Date$
% See also  iData/getaxis, iData/setaxis

a = circshift(a, varargin{:});
  
if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),a);
end
