function a = ylim(a, lims)
% b = ylim(s,[ ymin ymax ]) : Reduce iData Y axis limits
%
%   @iData/ylim function to reduce the Y axis (rank 2, rows) limits
%     ylim(s) returns the current Y axis limits. 
%
% input:  s: object or array (iData)
%         limits: new axis limits (vector)
% output: b: object or array (iData)
% ex:     b=ylim(a);
%
% See also iData, iData/plot, iData/ylabel

axisvalues = getaxis(a, 2);
if isempty(axisvalues), return; end
if nargin == 1
  a=[ min(axisvalues) max(axisvalues) ]; 
  return
end

index=find(lims(1) <= axisvalues & axisvalues <= lims(2));
s.type='()';
s.subs={ ':', index };
a = subsref(a,s);

