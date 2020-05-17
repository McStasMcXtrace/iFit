function a = zlim(a, lims, exclude)
% ZLIM Z limits (3rd axis).
%   ZL = ZLIM(A)            gets the Z limits.
%   Undefined axis returns [NaN NaN] as limits.
%
%   ZLIM(A,[ZMIN ZMAX])     sets the Z limits. 
%
%   ZLIM(A,[ZMIN ZMAX], 'exclude') removes the specified range instead of keeping it.
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plot, estruct/zlabel

% handle input estruct arrays
if nargin < 2, lims = ''; end
if nargin < 3, exclude = ''; end
if numel(a) > 1
  s = [];
  for index=1:numel(a)
    s = [ s ; feval(mfilename, a(index), lims, exclude) ];
  end
  if ~isempty(lims)
    a = reshape(s, size(a));
  else
    a = s;
  end
  return
end

axisvalues = getaxis(a, 3);
if isempty(axisvalues), a=[nan nan]; return; end
if isempty(lims)
  a=[ min(axisvalues(:)) max(axisvalues(:)) ]; 
  return
end

if ~isempty(exclude)
  index=find(lims(1) > axisvalues | axisvalues > lims(2));
else
  index=find(lims(1) < axisvalues & axisvalues < lims(2));
end
s.type='()';
s.subs={ ':', ':', index };
cmd=a.Command;
a = subsref(a,s);
a.Command=cmd;
a=history(a, mfilename, a, lims);

