function result = Math_1_unary
  a = iData([ ifitpath 'Data/ILL_IN6.dat' ]);
  op = {'asin', 'acos','atan','cos','sin','exp','log','log10','sqrt','tan','transpose',...
    'ctranspose', 'sparse','full', 'floor','ceil','round',...
    'del2','asinh','atanh','acosh','sinh','cosh','tanh', ...
    'sign','isfinite','isnan','isinf',...
    'isscalar','isvector','issparse','isreal','isfloat','isnumeric','isinteger','islogical',...
    'uminus','abs','real','imag','uplus','not',...
    'flipud','fliplr','permute','conj','cumtrapz','sum','prod','trapz','cumsum','cumprod'};
  result = [ 'OK  log floor sqrt cos ... (' num2str(length(op)) ' operators)' ];
  failed = '';
  for index=1:length(op)
    try
      feval(op{index}, a);
    catch
      failed = [ failed ' ' op{index} ];
    end
  end
  if length(failed)
    result = [ 'FAILED ' failed ];
  end
