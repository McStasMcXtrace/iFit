function a=shiftdim(a, varargin)
% SHIFTDIM Shifts an object axis of specified rank by a value.
%   A=SHIFTDIM(A, RANK, SHIFT) shifts object axes RANK by SHIFT.
%   This is equivalent to A{RANK} = A{RANK}+SHIFT; and MOVEAXIS.
%
% Example: a=iData(peaks); shiftdim(a, 1, 5); max(getaxis(a,1))==54
% Version: $Date$ $Version$ $Author$
% See also  iData/getaxis, iData/setaxis

a = moveaxis(a, varargin{:});
