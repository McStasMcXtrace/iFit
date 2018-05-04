function a = zlim(a, lims, exclude)
% b = zlim(s,[ zmin zmax ]) : Reduce iData Z axis limits
%
%   @iData/zlim function to reduce the Z axis (rank 3) limits
%     zlim(s) returns the current Z axis limits. 
%     Undefined axis returns [NaN NaN] as limits.
%
%   zlim(s, [min max], 'exclude') removes the specified range instead of keeping it.
%
% input:  s: object or array (iData)
%         limits: new axis limits (vector)
% output: b: object or array (iData)
% ex:     b=zlim(a);
%
% Version: $Date$
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
a=iData_private_history(a, mfilename, a, lims);

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),a);
end
