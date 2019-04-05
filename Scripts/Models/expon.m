function y=expon(varargin)
% y = expon(p, x, [y]) : Exponential decay
%
%   iFunc/expon Exponential decay fitting function
%     p(2)=Tau is the exponential decay parameter, in inverse 'x' units.
%     y=p(3)+p(1)*exp(-x/p(2));
%
% expon(decay)          creates a model with specified decay constant
% expon([ parameters ]) creates a model with specified model parameters
%
% input:  p: Exponential decay model parameters (double)
%            p = [ Amplitude Tau BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=expon([1 0.5 0], -0:10); or plot(expon)
%
% Version: $Date$ $Version$ $Author$
% See also iData, iFunc/fits, iFunc/plot, twoexp, sinedamp
% 

y.Name           = [ 'Exponential decay (1D) [' mfilename ']' ];
y.Description    = 'Exponential decay';
y.Parameters     = {'Amplitude','Tau decay in inverse "x" unit', 'Background'};
y.Expression     = @(p,x) p(3)+p(1)*exp(-x/p(2));
y.Dimension      = 1;         % dimensionality of input space (axes) and result

% use ifthenelse anonymous function
% <https://blogs.mathworks.com/loren/2013/01/10/introduction-to-functional-programming-with-anonymous-functions-part-1/>
% iif( cond1, exec1, cond2, exec2, ...)
iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
y.Guess     = @(x,y) iif(...
  ~isempty(y), @() [ ...
   exp(subsref(polyfit(x(:),log(y(:)-min(y(:))+0.01*abs(min(y(:)))),1), struct('type','()', 'subs',{{2}}))) ...
    -1/(subsref(polyfit(x(:),log(y(:)-min(y(:))+0.01*abs(min(y(:)))),1), struct('type','()', 'subs',{{1}}))-(abs(subsref(polyfit(x(:),log(y(:)-min(y(:))+0.01*abs(min(y(:)))),1), struct('type','()', 'subs',{{1}}))) < 1e-2)*.1) ...
    min(y(:)) ], ...
  true            , @() [1 0.5 0]);
y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 varargin{1} 0 ]};
  end
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end


