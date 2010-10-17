function a = xlabel(a, lab)
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
% Version: $Revision: 1.4 $
% See also iData, iData/plot, iData/label, iData/ylabel, iData/zlabel, iData/clabel

if nargin ==1
	a = label(a, 1);
	return
else
	a = label(a, 1, lab);
end

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),a);
end

