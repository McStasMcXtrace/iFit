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
    y=p(4)+p(1)*exp(-0.5*((x-p(2))/p(3)).^2);
else
  if nargin==0, action='identify'; end
  if ischar(p),     action=p; 
  elseif ischar(x), action=x; 
  elseif isstruct(p) & isfield(p,'action'), action = p.action; end
  
  switch action
  case 'identify'
    y.Name           = 'Gaussian (1D)';
    y.Parameters     = {'Amplitude','Centre','HalfWidth','Background'};
    y.Dimension      = 1;         % dimensionality of input space (axes)
    y.Guess          = [1 0 1 1]; % default parameters
  end
end

