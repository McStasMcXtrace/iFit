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
% Version: $Revision: 1.17 $
% See also iData, iData/plus, iData/prod, iData/cumcat, iData/mean
if nargin == 1 & isa(dim, 'iData') & length(dim) >= 1 % syntax: cat([a])
  s = cat(1, dim);
  return
end

if ~isa(a, 'iData')
  iData_private_error(mfilename,['syntax is cat(dim, iData, ...)']);
end

if length(varargin) > 1  % syntax: cat(dim,a,b,c,...)
  if numel(a) == 1, s=a; 
  else s=a(:); end
  for index=1:length(varargin)
    s = [ s ; varargin{index} ];
  end
  clear varargin
  s = cat(dim, s);
  return
end
% syntax is now: cat(dim,[a(:)])
a=a(:);
if length(a) <= 1, s=a; return; end
if dim <= 0, dim=1; end

% removes warnings during interp
iData_private_warning('enter', mfilename);

% syntax is now: cat(dim,[a b c ... ])
% first need to compute union axes, but not for dimension 'dim'
c_axis = iData_private_caxis(a,'union');
% use common axes on all axes except dim
for index=1:numel(a)
  x = getaxis(a(index), dim);
	if length(x) == 1 || length(x) == length(a)
		c_axis{dim} = x; % restore initial 'dim' axis from object. Others are the common axes.
	end
	a(index) = interp(a(index), c_axis);
end

% now catenate Signal, Error and Monitor 
lab = label(a, 0); if iscellstr(lab), lab=[ lab{1} '...' ]; end
for index=1:numel(a)
  s{index}=get(a(index),'Signal');
  if isvector(s{index}), ss=s{index}; ss=ss(:); s{index}=ss; end
end
ss = cat(dim, s{:});

% Error handling
% first test if all Errors are sqrt(this.Signal) (that is default=[])
% first test if all Monitors are 1 (that is default=[])
se = 1; sm = 1;
for index=1:numel(a)
  if ~strcmp(getalias(a(index),'Error'), 'sqrt(this.Signal)')
    se = 0;
  end
  if ~isnumeric(getalias(a(index), 'Monitor')) || length(getalias(a(index), 'Monitor')) > 1
    sm = 0;
  end
end

% then decide what to catenate for Error
if se == 1  % all Errors are default, just copy the default
  se = [];
else
  % some Errors are not default: we catenate all of them as values
  for index=1:length(a)
    se = getalias(a(index),'Error');
    if ~isnumeric(se) || ~isscalar(se)
      se = get(a(index),'Error');
    end
    s{index} = ones(size(get(a(index),'Signal'))).*se;
    if isvector(s{index}), se=s{index}; se=se(:); s{index}=se; end
  end
  se = cat(dim, s{:});
end

% then decide what to catenate for Monitors
if sm == 1  % all Monitors are default, just copy the default
  sm = [];
else
  % some Monitors are not default: we catenate all of them as values
  for index=1:numel(a)
    sm = getalias(a(index),'Monitor');
    if ~isnumeric(sm) || ~isscalar(sm)
      sm = get(a(index),'Monitor');
    end
    s{index} = ones(size(get(a(index),'Signal'))).*sm;
    if isvector(s{index}), sm=s{index}; sm=sm(:); s{index}=sm; end
  end
  sm = cat(dim, s{:});
end

% and extend axis 'dim'
for index=1:numel(a)
  s{index}=getaxis(a(index),dim); 
  if isempty(s{index}), s{index}=index;
  elseif isvector(s{index}), sx=s{index}; sx=sx(:); s{index}=sx; end
end

if isvector(s{1})
  sx = cat(1, s{:});
else
  sx = cat(dim, s{:});
end

% now build final result
cmd = get(a(1),'Command');
s = copyobj(a(1));  % with extended (union) axes
setalias(s,'Signal', ss, [ 'catenated ' lab ]);     % Store Signal
setalias(s,'Error',  se);
setalias(s,'Monitor',sm);
% create or modify catenated dimension axis
dx=getaxis(a(1),num2str(dim));
if isempty(dx)
  dx=[ 'Axis_' num2str(dim) ];
end
setaxis(s, dim, dx, sx);
s.Command=cmd;
s = iData_private_history(s, mfilename, dim, a(1), a(2));

% reset warnings during interp
iData_private_warning('exit', mfilename);


