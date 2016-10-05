function p = iFunc_feval_guess(this, varargin)
% private function to evaluate a guess in a reduced environment so that 
% internal function variables do not affect the result. 
% Guess=char as fhandle are handled directly in the core function
  ax = 'x y z t u v w';
  p  = [];
  this.Dimension = abs(this.Dimension);
  
  if this.Dimension
    eval([ '[' ax(1:(2*this.Dimension)) ']=deal(varargin{' mat2str(1:this.Dimension) '});' ]);
  end
  if length(varargin) > this.Dimension && ~isempty(varargin{this.Dimension+1}) && isnumeric(varargin{this.Dimension+1})
    signal = varargin{this.Dimension+1};
  else
    signal = 1;
  end
  clear ax
  % moments of distributions (used in some Guesses, e.g. gauss, lorz, ...)
  m1 = @(x,s) sum(s(:).*x(:))/sum(s(:));
  m2 = @(x,s) sqrt(abs( sum(x(:).*x(:).*s(:))/sum(s(:)) - m1(x,s).^2 ));
  try
    p = eval(this.Guess);     % returns model vector with output arg
  end
  if isempty(p)
    try
      eval(this.Guess);       % no output arg: returns model vector and redefines 'p'
    catch
      p = [];
    end
  end
