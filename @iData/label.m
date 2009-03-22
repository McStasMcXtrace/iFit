function a = label(a, alias, lab)
% b = label(s, alias, label) : Change iData label for a given alias/axis
%
%   @iData/label function to set/get labels
%     label(s, alias) returns the current label
%   The input iData object is updated if no output argument is specified.
%
% input:  s: object or array (iData)
%         alias: name of the alias or index of the axis (char/numeric)
%         label: new label (char/cellstr)
% output: b: object or array (iData)
% ex:     b=label(a,'x','new xlabel'); b=label(a,'x'); b=label(a, 1,'new xlabel');
%
% Version: $Revision: 1.2 $
% See also iData, iData/plot, iData/xlabel, iData/ylabel, iData/zlabel

if nargin < 2, alias=[]; end
if isempty(alias), return; end
if isnumeric(alias)
	[link, lab0] = getaxis(a, num2str(alias));
else
	[link, lab0] = getalias(a, alias);
end

if nargin == 2
  a=lab0; 
  return
end
cmd=a.Command;
setalias(a, link, getalias(a, link), lab);
a.Command=cmd;
a = iData_private_history(a, mfilename, a, alias, lab);

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),a);
end

