function a = xlabel(a, label)
% b = xlabel(s,label) : Change iData X axis label
%
%   @iData/xlabel function to change the X axis (rank 1, columns) label
%     xlabel(s) returns the current X axis label
%   The input iData object is updated if no output argument is specified.
%
% input:  s: object or array (iData)
%         label: new X label (char/cellstr)
% output: b: object or array (iData)
% ex:     b=xlabel(a);
%
% See also iData, iData/plot

[axisdef, lab] = getaxis(a, '1');
if isempty(axisdef), return; end
if nargin == 1
  a=lab; 
  return
end

setalias(a, axisdef, getalias(a, axisdef), label);

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),a);
end

