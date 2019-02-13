function sab=Sqw2Sab(self)
% convert a Sqw2D into a Sab object

% the final object will have alpha beta axes. These axes must be changed into qw
% then evaluate the original model, and set back axes.

% new parameters:
%  Temperature
%  Weight
%  

  sab=copyobj(self);
  index=numel(sab.Parameters)+1;
  
   if isvector(sab.Guess) && isnumeric(sab.Guess)
    sab.Guess(index)   = 300;
    sab.Guess(index+1) = 12;
  else
    if ~iscell(sab.Guess), sab.Guess = { sab.Guess }; end
    sab.Guess{end+1} = [ 300 12 ];
  end
  
  sab.Parameters{end+1} = 'Temperature [K]'; 
  sab.Parameters{end+1} = 'Mass Material molar weight [g/mol]';
  sab.Expression{end+1} = [ 'T = p(' num2str(index) '); M=p(' num2str(index+1) ');' ];
  
  sab=iFunc(sab);
  
  % build the expression
  sab.Expression{end+1} = [ 'alpha = x; beta = y;' ];
  sab.Expression{end+1} = 'q2toE   = 2.072/M; % [Angs^-2] to [meV] ';
  % the definition of alpha states that it is per M(atom/molecule), so we divide
  sab.Expression{end+1} = 'C       = 11.6003/T;';
  sab.Expression{end+1} = 'q =  sqrt(alpha/(C*q2toE/M));';
  % sab.Expression{end+1} = 'w = -beta/C;'; % not needed
  % Jacobian for S(q,w) -> S(a,b) is J=(dq.dw)/(dalpha.dbeta) = 1/(2*q*C^2*q2toE)
  sab.Expression{end+1} = 'signal = signal.*q*(2*C*C*q2toE);';
  sab.Expression{end+1} = 'x = alpha; y = beta;';
  sab.Expression{end+1} = 'signal(imag(signal) | imag(q))=0;';
  sab.Expression{end+1} = 'signal=real(signal);';
  
  sab=iFunc(sab);
