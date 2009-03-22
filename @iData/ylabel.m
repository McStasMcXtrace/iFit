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
% Version: $Revision: 1.3 $
% See also iData, iData/plot, iData/xlabel, iData/zlabel, iData/clabel

if nargin ==1
	a = label(a, 2);
	return
else
	a = label(a, 2, lab);
end

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),a);
end

