function y=expstretched(varargin)
% y = expstretched(p, x, [y]) : Exponential decay
%
%   iFunc/expstretched Stretched exponential decay fitting function
%     Tau is the expeonential decay parameter, in inverse 'x' units.
%     y=p(4)+p(1)*exp(-(x/p(2)).^p(3));
%
% expstretched(decay)          creates a model with specified decay constant
% expstretched([ parameters ]) creates a model with specified model parameters
%
%   This model can be used as a Debye-Waller factor, e.g. exp(-u.^2 q.^2/3)
%   <https://en.wikipedia.org/wiki/Debye%E2%80%93Waller_factor>
%   then fix p(1)=1; p(3)=2; p(4)=0; 
%   and  use p(2)=sqrt(3/<u^2>) where <u^2> is the mean squared displacement.
%
% input:  p: Stretched exponential decay model parameters (double)
%            p = [ Amplitude Tau Exponent BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=expstretched([1 0.2 2 0], -10:10); or plot(expstretched)
%
%         I will not buy this exponential; it is stretched.
%         <http://en.wikipedia.org/wiki/Dirty_Hungarian_Phrasebook>
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot
% (c) E.Farhi, ILL. License: EUPL.

y.Name           = [ 'Stretched Exponential decay (1D) [' mfilename ']' ];
y.Description    = 'Stretched Exponential decay';
y.Parameters     = {'Amplitude','Tau decay in inverse "x" unit', 'Exponent', 'Background'};
y.Expression     = @(p,x) p(4)+p(1)*exp(-(x/p(2)).^p(3));
y.Dimension      = 1;         % dimensionality of input space (axes) and result

% use ifthenelse anonymous function
% <https://blogs.mathworks.com/loren/2013/01/10/introduction-to-functional-programming-with-anonymous-functions-part-1/>
% iif( cond1, exec1, cond2, exec2, ...)
iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
y.Guess     = @(x,y) iif(...
  ~isempty(y), @() [ ...
   exp(subsref(polyfit(x(:),log(y(:)-min(y(:))+0.01*abs(min(y(:)))),1), struct('type','()', 'subs',{{2}}))) ...
    -1/(subsref(polyfit(x(:),log(y(:)-min(y(:))+0.01*abs(min(y(:)))),1), struct('type','()', 'subs',{{1}}))- ...
       (abs(subsref(polyfit(x(:),log(y(:)-min(y(:))+0.01*abs(min(y(:)))),1), struct('type','()', 'subs',{{1}}))) < 1e-2)*.1) ...
    1 min(y(:)) ], ...
  true            , @() [1 0.2 2 0]);

y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 varargin{1} 1 0 ]};
  end
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end

