function a = xlim(a, lims, exclude)
% XLIM X limits (2nd axis for ndims>=2).
%   XL = XLIM(A)            gets the X limits.
%   Undefined axis returns [NaN NaN] as limits.
%
%   XLIM(A,[XMIN XMAX])     sets the X limits. 
%
%   XLIM(A,[XMIN XMAX], 'exclude') removes the specified range instead of keeping it.
%
% Example: a=iData(peaks); b=xlim(a,[5 35]); all(xlim(b)==[5 35])
% Version: $Date$ $Version$ $Author$
% See also iData, iData/plot, iData/xlabel

% handle input iData arrays
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

if ndims(a) == 1
  axisvalues = getaxis(a, 1);
else
  axisvalues = getaxis(a, 2);
end
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
if ndims(a) > 1
  s.subs={ ':', index };
else
  s.subs={ index };
end
cmd=a.Command;
a = subsref(a,s);
a.Command=cmd;
a=history(a, mfilename, a, lims);
