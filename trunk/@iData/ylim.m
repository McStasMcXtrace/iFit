function a = ylim(a, lims)
% b = ylim(s,[ ymin ymax ]) : Reduce iData Y axis limits
%
%   @iData/ylim function to reduce the Y axis (rank 1, rows) limits
%     ylim(s) returns the current Y axis limits. 
%
% input:  s: object or array (iData)
%         limits: new axis limits (vector)
% output: b: object or array (iData)
% ex:     b=ylim(a);
%
% Version: $Revision: 1.5 $
% See also iData, iData/plot, iData/ylabel

axisvalues = getaxis(a, 1);
if isempty(axisvalues), return; end
if nargin == 1
  a=[ min(axisvalues) max(axisvalues) ]; 
  return
end

index=find(lims(1) <= axisvalues & axisvalues <= lims(2));
s.type='()';
if ndims(a) > 1
  s.subs={ index, ':' };
else
  s.subs={ index };
end
cmd=a.Command;
a = subsref(a,s);
a.Command=cmd;
a=iData_private_history(a, mfilename, a, lims);

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),a);
end
