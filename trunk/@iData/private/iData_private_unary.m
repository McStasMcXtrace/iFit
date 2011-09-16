function b = iData_private_unary(a, op)
% iData_private_unary: handles unary operations
%
% 'isscalar','isvector','issparse','isreal','isfloat','isnumeric','isinteger','islogical'
% 'asin', 'acos','atan','cos',sin','exp','log','log10','sqrt','tan','transpose'
% 'ctranspose', 'sparse','full', 'floor','ceil','round','sign','isfinite','isnan','isinf'
% 'del2'
% 'sign','isfinite','isnan','isinf'
% 'isscalar','isvector','issparse','isreal','isfloat','isnumeric','isinteger','islogical'
% 'uminus','abs','real','imag','uplus','not'
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
cmd=a.Command;
b = copyobj(a);
s = get(b,'Signal');
[dummy, sl] = getaxis(b, '0');

% operation on Error
e = get(b,'Error');
m = get(b,'Monitor');

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
new_s = feval(op, s);

switch op
case 'acos'
	e = -e./sqrt(1-s*s);
case 'asin'
	e = e./sqrt(1-s*s);
case 'atan'
	e = e./(1+s*s);
case 'cos'
	e = -e.*sin(s);
case 'exp'
	e = e.*exp(s);
case 'log'
	e = e./s;
case 'log10'
	e = e./(log(10)*s);
case 'sin'
	e = e.*cos(s);
case 'sqrt'
	e = e/(2*sqrt(s));
  m = m.^0.5;
case 'tan'
	c = cos(s);
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
case {'sparse','full'}
  % apply same operator on error and Monitor
	e = feval(op, e);
	m = feval(op, m);
case {'floor','ceil','round'}	
	% apply same operator on error
	e = feval(op, e);
case 'del2'
  s = 2*ndims(a);
  e = feval(op, e)*2*ndims(a);
case {'sign','isfinite','isnan','isinf'}
	% error should become zero (logical)
	e = zeros(size(s));
case {'isscalar','isvector','issparse','isreal','isfloat','isnumeric','isinteger','islogical','double','single','logical','find'}
	% result is a single value
	b = new_s;
	return
case {'uminus','abs','real','imag','uplus','not'}
	% retain error, do nothing
otherwise
  iData_private_error('unary',['Can not apply operation ' op ' on object ' a.Tag ]);
end

% update object
b = set(b, 'Signal', new_s, 'Error', abs(e), 'Monitor', m);
% test if we could update signal as expected, else we store the new value directly in the field
if ~isequal(get(b,'Signal'), new_s)
  b = setalias(b, 'Signal', new_s, [  op '(' sl ')' ]);
end
if ~isequal(get(b,'Error'), e)
  b = setalias(b, 'Error', e);
end
if ~isequal(get(b,'Monitor'), m)
  b = setalias(b, 'Monitor', m);
end
b.Command=cmd;
b = iData_private_history(b, op, a);  

% other methods to treat specifically
% diff, min, max, sum, prod, sort, trapz

