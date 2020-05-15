function a=moveaxis(a, dim, shift)
% MOVEAXIS Shifts an object axis of specified rank by a value.
%   A=MOVEAXIS(A, RANK, SHIFT) shifts object axes RANK by SHIFT.
%   This is equivalent to A{RANK} = A{RANK}+SHIFT;
%
% Example: a=estruct(peaks); moveaxis(a, 1, 5); max(getaxis(a,1))==54
% Version: $Date$ $Version$ $Author$
% See also  estruct/getaxis, estruct/setaxis

if nargin < 3
  return
end

% handle input estruct arrays
if numel(a) > 1
  for index=1:numel(a)
    a(index) = feval(mfilename, a(index), dim, shift);
  end
  return
end

a = setaxis(a, dim, getaxis(a, dim)+shift);
