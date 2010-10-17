function s_out = rmaxis(a_in,indexes)
% s = rmaxis(s, Axis) : delete iData axes
%
%   @iData/rmaxis function to delete iData axes.
%     As Axes are often required for the Signal to be used, 
%     default axes will be automatically re-created when needed.
%   The function works also when Axis is given as a cell.
%   The Axis can be specified as a rank (1-ndims) or an axis name. 
%   When the Axis is empty, all axes are removed.
%   The input iData object is updated if no output argument is specified.
%   Axis 1 is often labeled as 'x' (on columns), 2 as 'y' (on rows), etc...
%
% input:  s: object or array (iData)
%         AxisIndex: rank or alias name of the axis
% output: s: array (iData)
% ex:     rmaxis(iData, 1) removes the 'x' axis (rank 1)
%
% Version: $Revision: 1.1 $
% See also iData, iData/getaxis, iData/get, iData/set, iData/setaxis

% EF 27/07/00 creation
% EF 23/09/07 iData implementation

if isnumeric(indexes)
  s_out=setaxis(a_in, indexes, []);
else
  s_out=setaxis(a_in, [], indexes);
end

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),s_out);
end
