function a = zlabel(a, lab)
% ZLABEL Z-axis object label (rank 3, pages).
%   ZLABEL(s, 'label') changes the Z axis (rank 3, pages) label
%
%   ZLABEL(s) returns the current Z axis label
%
% Example: s=iData(1:10); zlabel(s,'Z'); strcmp(zlabel(s),'Z')
% Version: $Date$ $Version$ $Author$
% See also iData, iData/plot, iData/xlabel, iData/ylabel, iData/label, iData/clabel

if nargin ==1
	a = label(a, 3);
	return
else
	a = label(a, 3, lab);
end
