function a = ylim(a, lims, exclude)
% YLIM Y limits (1st axis for ndims>=2).
%   YL = YLIM(A)            gets the Y limits.
%   Undefined axis returns [NaN NaN] as limits.
%
%   YLIM(A,[YMIN YMAX])     sets the Y limits. 
%
%   YLIM(A,[YMIN YMAX], 'exclude') removes the specified range instead of keeping it.
%
% Example: a=iData(peaks); b=ylim(a,[5 35]); all(ylim(b)==[5 35])
% Version: $Date$ $Version$ $Author$
% See also iData, iData/plot, iData/ylabel

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
  if nargout == 0 & nargin == 2 & length(inputname(1))
    assignin('caller',inputname(1),a);
  end
  return
end

axisvalues = getaxis(a, 1);
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
if ndims(a) > 1 && numel(axisvalues) == max(size(axisvalues))
  s.subs={ index, ':' };
else
  s.subs={ index }; % this creates an 'event' object
end

cmd=a.Command;
a = subsref(a,s);
a.Command=cmd;
a=history(a, mfilename, a, lims);

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),a);
end
