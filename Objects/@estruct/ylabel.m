function l = ylabel(a, lab)
% YLABEL Y-axis object label (rank 1, rows).
%   YLABEL(s, 'label') changes the Y axis (rank 1, rows) label
%
%   YLABEL(s) returns the current Y axis label
%
% Example: s=estruct(1:10); ylabel(s,'Y'); strcmp(ylabel(s),'Y')
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plot, estruct/xlabel, estruct/label, estruct/zlabel, estruct/clabel

if nargin ==1
  if isvector(a) == 1
    l = label(a, 0);
  else
	  l = label(a, 1);
  end
  return
else
	if isvector(a) == 1
    l = label(a, 0, lab);
  else
	  l = label(a, 1, lab);
  end
end

