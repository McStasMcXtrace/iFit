function a = zlim(a, lims)
% b = zlim(s,[ zmin zmax ]) : Reduce iData Z axis limits
%
%   @iData/zlim function to reduce the Z axis (rank 3) limits
%     zlim(s) returns the current Z axis limits. 
%
% input:  s: object or array (iData)
%         limits: new axis limits (vector)
% output: b: object or array (iData)
% ex:     b=zlim(a);
%
% Version: $Revision: 1.2 $
% See also iData, iData/plot, iData/ylabel

axisvalues = getaxis(a, 3);
if isempty(axisvalues), return; end
if nargin == 1
  a=[ min(axisvalues) max(axisvalues) ]; 
  return
end

index=find(lims(1) <= axisvalues & axisvalues <= lims(2));
s.type='()';
s.subs={ ':', ':', index };
a = subsref(a,s);

