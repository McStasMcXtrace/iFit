function a = clim(a, lims, exclude)
% CLIM C limits (4th axis).
%   CL = CLIM(A)             gets the C limits.
%   Undefined axis returns [NaN NaN] as limits.
%
%   CLIM(A,[CMIN CMAX])     sets the C limits. 
%
%   CLIM(A,[CMIN CMAX], 'exclude') removes the specified range instead of keeping it.
%
% Example: a=iData(flow); all(isnan(clim(a)))
% Version: $Date$ $Version$ $Author$
% See also iData, iData/plot, iData/clabel

% handle input arrays
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

axisvalues = getaxis(a, 4);
if isempty(axisvalues), a=[nan nan]; return; end
if isempty(lims)
  a=[ min(axisvalues(:)) max(axisvalues(:)) ];
  return
end

if ~isempty(exclude)
  index=find(lims(1) >= axisvalues | axisvalues >= lims(2));
else
  index=find(lims(1) <= axisvalues & axisvalues <= lims(2));
end
s.type='()';
s.subs={ ':', ':', ':', index };
cmd=a.Command;
a = subsref(a,s);
a.Command=cmd;
a=history(a, mfilename, a, lims);

