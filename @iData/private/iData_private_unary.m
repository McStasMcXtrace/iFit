function b = iData_private_unary(a, op)
% iData_private_unary: handles unary operations

% handle input iData arrays
if length(a(:)) > 1
  b = a;
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

% operation on signal
s = feval(op, s);

% operation on Error
e = get(b,'Error');
m = get(b,'Monitor');
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
    m = m.^2;
case 'tan'
	c = cos(s);
	e = e./(c.*c);
case { 'transpose', 'ctranspose'}; % .' and ' respectively
	e = e';
	if ndims(b) > 1
  	x1 = getaxis(b, '1'); % axis names
  	x2 = getaxis(b, '2');
  	v1 = getaxis(b, 1);   % axis values
  	v2 = getaxis(b, 2);
  	if ~isempty(x2), b= setaxis(b, 1, x2, transpose(v2)); end
  	if ~isempty(x1), b= setaxis(b, 2, x1, transpose(v1)); end
  end
case {'sparse','full'}
  % apply same operator on error
	e = feval(op, e);
	b = set(b, 'Monitor', feval(op, get(b,'Monitor')));
case {'floor','ceil','round'}	
	% apply same operator on error
	e = feval(op, e);
case {'sign','isfinite','isnan','isinf'}
	% error should become zero (logical)
	e = zeros(size(s));
case {'isscalar','isvector','issparse','isreal','isfloat','isnumeric','isinteger','islogical','double','single','logical','find'}
	% result is a single value
	b = s;
	return
case {'uminus','abs','real','imag','uplus','not'}
	% retain error, do nothing
otherwise
  iData_private_error('unary',['Can not apply operation ' op ' on object ' a.Tag ]);
end

% update object
b = set(b, 'Signal', s, 'Error', abs(e));
b = setalias(b, 'Signal', s, [  op '(' sl ')' ]);
b.Command=cmd;
b = iData_private_history(b, op, b);  

% other methods to treat specifically
% diff, min, max, sum, prod, sort, trapz

