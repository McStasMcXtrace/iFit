function y=gauss(p, x)
% y = gauss(p, x) : Gaussian
%
%   iFunc/gauss Gaussian fitting function
%   the function called with a char argument performs specific actions.
%
% input:  p: Gaussian model parameters
%            p = [ Amplitude Centre HalfWidth BackGround ]
%         x: axis (double) or action e.g. 'identify' (char)
% output: y: model value
% ex:     y=gauss([1 0 1 1], -10:10); or y=gauss('identify')

if nargin==2
    if isempty(p) 
      id=feval(mfilename, 'identify');
      p=id.Guess;
    end
    y=p(4)+p(1)*exp(-0.5*((x-p(2))/p(3)).^2);
else
  action='';
  if nargin == 1 & ischar(p), action=p; end
  if isempty(action), action='identify'; end
  
  switch action
  case 'identify'
    y.Type           = 'iFit fitting function';
    y.Name           = 'Gaussian (1D)';
    y.Parameters     = {'Amplitude','Centre','HalfWidth','Background'};
    y.Dimension      = 1;         % dimensionality of input space (axes) and result
    y.Guess          = [1 0 1 1]; % default parameters
  end
end

