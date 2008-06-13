function a = clabel(a, label)
% b = clabel(s,label) : Change iData C axis label
%
%   @iData/clabel function to change the C axis (rank 4) label
%     clabel(s) returns the current C axis label
%   The input iData object is updated if no output argument is specified.
%
% input:  s: object or array (iData)
%         label: new C label (char/cellstr)
% output: b: object or array (iData)
% ex:     b=clabel(a);
%
% Version: $Revision: 1.2 $
% See also iData, iData/plot

[axisdef, lab] = getaxis(a, '4');
if isempty(axisdef), return; end
if nargin == 1
  a=lab; 
  return
end

setalias(a, axisdef, getalias(a, axisdef), label);

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),a);
end

