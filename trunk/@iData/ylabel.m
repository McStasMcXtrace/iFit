function a = ylabel(a, label)
% b = ylabel(s,label) : Change iData Y axis label
%
%   @iData/ylabel function to change the Y axis (rank 2, rows) label
%     ylabel(s) returns the current Y axis label
%   The input iData object is updated if no output argument is specified.
%
% input:  s: object or array (iData)
%         label: new Y label (char/cellstr)
% output: b: object or array (iData)
% ex:     b=ylabel(a);
%
% See also iData, iData/plot

[axisdef, lab] = getaxis(a, '2');
if isempty(axisdef), return; end
if nargin == 1
  a=lab; 
  return
end

setalias(a, axisdef, getalias(a, axisdef), label);

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),a);
end

