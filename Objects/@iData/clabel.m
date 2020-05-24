function a = clabel(a, lab)
% CLABEL 4-th axis object labels ('Color' axis label)
%   CLABEL(s, 'label') changes the C axis (rank 4) label
%
%   CLABEL(s) returns the current C axis label
%
% Example: s=iData(1:10); clabel(s,'C'); strcmp(clabel(s),'C')
% Version: $Date$ $Version$ $Author$
% See also iData, iData/plot, iData/xlabel, iData/ylabel, iData/zlabel, iData/label

if nargin ==1
	a = label(a, 4);
	return
else
	a = label(a, 4, lab);
end
