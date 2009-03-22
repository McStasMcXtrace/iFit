function a = xlim(a, lims)
% b = xlim(s,[ xmin xmax ]) : Reduce iData X axis limits
%
%   @iData/xlim function to reduce the X axis (rank 1, columns) limits
%     xlim(s) returns the current X axis limits. 
%
% input:  s: object or array (iData)
%         limits: new axis limits (vector)
% output: b: object or array (iData)
% ex:     b=xlim(a);
%
% Version: $Revision: 1.3 $
% See also iData, iData/plot, iData/xlabel

axisvalues = getaxis(a, 1);
if isempty(axisvalues), return; end
if nargin == 1
  a=[ min(axisvalues) max(axisvalues) ]; 
  return
end

index=find(lims(1) < axisvalues & axisvalues < lims(2));
s.type='()';
if ndims(a) > 1
  s.subs={ index , ':' };
else
  s.subs={ index };
end
cmd=a.Command;
a = subsref(a,s);
a.Command=cmd;
a=iData_private_history(a, mfilename, a, lims);
