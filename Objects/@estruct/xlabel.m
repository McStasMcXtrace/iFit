function l = xlabel(a, lab)
% XLABEL X-axis object label (rank2, columns).
%   XLABEL(s, 'label') changes the X axis (rank 2, columns) label
%
%   XLABEL(s) returns the current X axis label
%
% Example: s=estruct(1:10); xlabel(s,'X'); strcmp(xlabel(s),'X')
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plot, estruct/label, estruct/ylabel, estruct/zlabel, estruct/clabel

if nargin ==1
  if isvector(a) == 1
	  l = label(a, 1);
  else
    l = label(a, 2);
  end
  return
else
  if isvector(a) == 1
    l = label(a, 1, lab);
  else
	  l = label(a, 2, lab);
  end
end
