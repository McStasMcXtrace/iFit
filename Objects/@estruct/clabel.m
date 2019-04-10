function a = clabel(a, lab)
% CLABEL 4-th axis object labels ('Color' axis label)
%   CLABEL(s, 'label') changes the C axis (rank 4) label
%
%   CLABEL(s) returns the current C axis label
%
% Example: s=estruct(1:10); clabel(s,'C'); strcmp(clabel(s),'C')
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plot, estruct/xlabel, estruct/ylabel, estruct/zlabel, estruct/label

if nargin ==1
	a = label(a, 4);
	return
else
	a = label(a, 4, lab);
end
