function b = iData_private_unary(a, op)
% iData_private_unary: handles unary operations
%
% 'asin', 'acos','atan','cos','sin','exp','log','log10','sqrt','tan','transpose'
% 'ctranspose', 'sparse','full', 'floor','ceil','round'
% 'asinh','atanh','acosh','sinh','cosh','tanh'
% 'del2'
% 'sign','isfinite','isnan','isinf'
% 'isscalar','isvector','issparse','isreal','isfloat','isnumeric','isinteger','islogical'
% 'uminus','abs','real','imag','uplus','not'
% 'flipud','fliplr'
%
% present but not used here: 'double','single','logical','find'

% handle input iData arrays
if length(a(:)) > 1
  switch op
  case {'isscalar','isvector','issparse','isreal','isfloat','isnumeric','isinteger','islogical'}
    b = ones(size(a));
  otherwise
    b =a;
  end
  for index=1:length(a(:))
    b(index) = iData_private_unary(a(index), op);
  end
  b = reshape(b, size(a));
  return
end

iData_private_warning('enter',[ mfilename ' ' op ]);

cmd=a.Command;
b = copyobj(a);
s = subsref(b,struct('type','.','subs','Signal'));
[dummy, sl] = getaxis(b, '0');

% operation on Error
e = subsref(b,struct('type','.','subs','Error'));
m = subsref(b,struct('type','.','subs','Monitor'));

% operation on signal
if strcmp(op, 'sparse')
  if ndims(a) > 2
    iData_private_error('unary',['Operation ' op ' can only be used on 2d data sets. Object ' a.Tag ' is ' num2str(ndims(a)) 'd.' ]);
  end
  if ~strcmp(class(s), 'double') && ~strcmp(class(s), 'logical')
    s = double(s);
  end
  if ~strcmp(class(e), 'double') && ~strcmp(class(e), 'logical')
    e = double(e);
  end
  if ~strcmp(class(m), 'double') && ~strcmp(class(m), 'logical')
    m = double(m);
  end
end

if ~isempty(find(strcmp(op, {'norm','asin', 'acos','atan','cos','sin','exp','log',...
 'log10','sqrt','tan','asinh','atanh','acosh','sinh','cosh','tanh'}))) ...
   && not(all(m(:) == 0 | m(:) == 1))
  s = genop(@rdivide, s, m);
  e = genop(@rdivide, e, m);
end

% non-linear operators should perform on the Signal/Monitor
% and then multiply again by the Monitor
new_s = feval(op, s); % new Signal value is set HERE <==========================

switch op
case 'acos'
	e = -e./sqrt(1-s.*s);
case 'acosh'
  e = e./sqrt(s.*s-1);
case 'asin'
	e = e./sqrt(1-s.*s);
case 'asinh'
  e = e./sqrt(1+s.*s);
case 'atan'
	e = e./(1+s.*s);
case 'atanh'
  e = e./(1-s.*s);
case 'cos'
	e = -e.*sin(s);
case 'cosh'
  e = e.*sinh(s);
case 'exp'
	e = e.*exp(s);
case 'log'
	e = e./s;
case 'log10'
	e = e./(log(10)*s);
case 'sin'
	e = e.*cos(s);
case 'sinh'
  e = e.*cosh(s);
case 'sqrt'
	e = e./(2*sqrt(s));
  m = m.^0.5;
case 'tan'
	c = cos(s);
	e = e./(c.*c);
case 'tanh'
  c = cosh(s);
  e = e./(c.*c);
case { 'transpose', 'ctranspose'}; % .' and ' respectively
	e = feval(op, e);
	m = feval(op, m);
	if ndims(b) > 1
  	x1 = getaxis(b, '1'); % axis names
  	x2 = getaxis(b, '2');
  	v1 = getaxis(b, 1);   % axis values
  	v2 = getaxis(b, 2);
  	if ~isempty(x2), b= setaxis(b, 1, x2, transpose(v2)); end
  	if ~isempty(x1), b= setaxis(b, 2, x1, transpose(v1)); end
  end
case {'sparse','full','flipud','fliplr'}
  % apply same operator on error and Monitor
	e = feval(op, e);
	m = feval(op, m);
case {'floor','ceil','round'}	
	% apply same operator on error
	e = feval(op, e);
case 'del2'
  new_s = new_s*2*ndims(a);
  e = 2*ndims(a)*del2(e);
case {'sign','isfinite','isnan','isinf'}
	b = new_s;
	return
case {'isscalar','isvector','issparse','isreal','isfloat','isnumeric','isinteger', ...
      'islogical','double','single','logical','find','norm'}
	% result is a single value
	b = new_s;
	iData_private_warning('exit',mfilename);
	return
case {'uminus','abs','real','imag','uplus','not'}
	% retain error, do nothing
otherwise
  iData_private_error('unary',['Can not apply operation ' op ' on object ' a.Tag ]);
end

if ~isempty(find(strcmp(op, {'norm','asin', 'acos','atan','cos','sin','exp','log',...
 'log10','sqrt','tan','asinh','atanh','acosh','sinh','cosh','tanh'}))) ...
   && not(all(m(:) == 0 | m(:) == 1))
  new_s = genop(@times, new_s, m);
  e     = genop(@times, e, m);
end

% update object
e = abs(e);
b = set(b, 'Signal', new_s, 'Error', e, 'Monitor', m);
% test if we could update signal as expected, else we store the new value directly in the field
if ~isequal(subsref(b,struct('type','.','subs','Signal')), new_s)
  b = setalias(b, 'Signal', new_s, [  op '(' sl ')' ]);
end
if ~isequal(subsref(b,struct('type','.','subs','Error')), e)
  b = setalias(b, 'Error', e);
end
if ~isequal(subsref(b,struct('type','.','subs','Monitor')), m)
  b = setalias(b, 'Monitor', m);
end
b.Command=cmd;
b = iData_private_history(b, op, a);  

% other methods to treat specifically
% diff, min, max, sum, prod, sort, trapz

iData_private_warning('exit',mfilename);

