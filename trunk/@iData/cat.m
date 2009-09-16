function s = cat(dim,a,varargin)
% s = cat(dim,a,...) : catenate iData objects elements along dimension
%
%   @iData/cat function to catenate iData objects elements along dimension dim
%     cat(dim,a,b,c) catenates along axis of rank dim. The axis is then extended.
%       cat(1,a) accumulates on first dimension (columns)
%
% input:  a: object or array (iData/array of)
%         dim: dimension to accumulate (int)
% output: s: catenated data set (iData)
% ex:     c=cat(1,a,b); c=cat(1,[ a b ]); 
%
% Version: $Revision: 1.7 $
% See also iData, iData/plus, iData/prod, iData/cumcat, iData/mean
if nargin == 1 & isa(dim, 'iData') & length(dim) >= 1 % syntax: cat([a])
  s = cat(1, dim);
  return
end

if ~isa(a, 'iData')
  iData_private_error(mfilename,['syntax is cat(dim, iData, ...)']);
end

if length(varargin) >= 1  % syntax: cat(dim,a,b,c,...)
  s=a(:);
  for index=1:length(varargin)
    s = [ s ; varargin{index} ];
  end
  s = cat(dim, s);
  return
end
% syntax is now: cat(dim,[a(:)])
a=a(:);
if length(a) <= 1, s=a; return; end

% removes warnings during interp
try
  warn.set = warning('off','iData:setaxis');
  warn.get = warning('off','iData:getaxis');
catch
  warn = warning('off');
end

% syntax is now: cat(dim,[a b c ... ])
% first need to compute union axes, but not for dimension 'dim'
c_axis = iData_private_caxis(a);

% use common axes on all axes except dim
for index=1:length(a)
	if length(getaxis(a(index), dim)) == 1 || length(getaxis(a(index), dim)) == length(a)
		c_axis{dim} = getaxis(a(index), dim); % restore initial 'dim' axis from object. Others are the common axes.
	end
	a(index) = interp(a(index), c_axis);
end

% now catenate Signal, Error and Monitor 
[link, label]          = getalias(a(1), 'Signal');
for index=1:length(a)
  s{index}=get(a(index),'Signal');
end
ss = cat(dim, s{:});
for index=1:length(a)
  s{index}=get(a(index),'Error');
end
se = cat(dim, s{:});
for index=1:length(a)
  s{index}=get(a(index),'Monitor');
end
sm = cat(dim, s{:});

% and extend axis 'dim'
for index=1:length(a)
  s{index}=getaxis(a(index),dim);
  if isempty(s{index}), s{index}=index; end
end
sx = cat(dim, s{:});

% now build final result
cmd = get(a(1),'Command');
s = copyobj(a(1));  % with extended (union) axes
setalias(s,'Signal', ss, [ 'catenated ' label ]);     % Store Signal
setalias(s,'Error',se);
setalias(s,'Monitor',sm);
dx=getaxis(a(1),num2str(dim));
if isempty(dx)
  dx=[ 'Axis_' num2str(dim) ];
end
setaxis(s, dim, dx, sx);
s.Command=cmd;
s = iData_private_history(s, mfilename, dim, a(1), a(2));

% reset warnings during interp
try
  warning(warn.set);
  warning(warn.get);
catch
  warning(warn);
end


