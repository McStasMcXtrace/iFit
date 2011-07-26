function a = xlim(a, lims)
% b = xlim(s,[ xmin xmax ]) : Reduce iData X axis limits
%
%   @iData/xlim function to reduce the X axis (rank 2, columns) limits
%     xlim(s) returns the current X axis limits. For 1D objects, the 
%
% input:  s: object or array (iData)
%         limits: new axis limits (vector)
% output: b: object or array (iData)
% ex:     b=xlim(a);
%
% Version: $Revision: 1.5 $
% See also iData, iData/plot, iData/xlabel

if ndims(a) == 1
  axisvalues = getaxis(a, 1);
else
  axisvalues = getaxis(a, 2);
end
if isempty(axisvalues), return; end
if nargin == 1
  a=[ min(axisvalues) max(axisvalues) ]; 
  return
end

index=find(lims(1) < axisvalues & axisvalues < lims(2));
s.type='()';
if ndims(a) > 1
  s.subs={ ':', index };
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
