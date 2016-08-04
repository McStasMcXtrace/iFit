function y=voigt(varargin)
% y = voigt(p, x, [y]) : Voigt function
%
%   iFunc/voigt Voigt fitting function, including Bose factor.
%     The Width parameters of the Gaussian and Lorentzian are the Half Widths.
%
% voigt(width)          creates a model with a specified width
% voigt([ parameters ]) creates a model with specified model parameters
%
% Reference: http://en.wikipedia.org/wiki/Voigt_profile
%
% input:  p: Voigt model parameters (double)
%            p = [ Amplitude Centre Width_Gauss Width_Lorz BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=voigt([1 0 1 1], -10:10); or plot(voigt);
%
% Version: $Date$
% See also iFunc, iFunc/fits, iFunc/plot
% (c) E.Farhi, ILL. License: EUPL.

y.Name      = [ 'Voigt (1D) [' mfilename ']' ];
y.Parameters={'Amplitude','Centre','Width_Gauss Half width',' Width_Lorz Half width', 'Background'};
y.Description='Single 1D Voigt model (gaussian x lorentzian)';
y.Dimension =1;
% MZ <mzinkin@sghms.ac.uk> adapted from DFM
y.Expression= { ...
  'p(3:4) = p(3:4) * 2;', ...
  'Nv = 16;', ...
	'Mv=2*Nv; M2=2*Mv;', ...
	'Lv=sqrt(Nv/sqrt(2));', ...
	'tt=(Lv*tan([-Mv+1:1:Mv-1]''*pi/M2)).^2;', ...
	'av=real(fft(fftshift([0; exp(-tt).*(Lv^2+tt)])))/M2;', ...
	'zz= -sqrt(log(2))/p(3)*p(4)*(1+2*i*(x-p(2)));', ...
	'lv=Lv-zz;', ...
	'signal=p(5)+p(1)*real(2*polyval(flipud(av(2:Nv+1)),(Lv+zz)./lv) ./lv.^2+(1/sqrt(pi))*ones(size(zz)) ./lv);'};

% moments of distributions
m1 = @(x,s) sum(s(:).*x(:))/sum(s(:));
m2 = @(x,s) sqrt(abs( sum(x(:).*x(:).*s(:))/sum(s(:)) - m1(x,s).^2 ));

y.Guess     = @(x,s) [ NaN m1(x, s-min(s(:))) m2(x, s-min(s(:))) m2(x, s-min(s(:))) NaN ];


y = iFunc(y);

if nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1
    varargin = {[ 1 0 varargin{1} varargin{1} 0]};
  elseif length(varargin{1}) == 2
    varargin = {[ 1 0 varargin{1} 0]};
  end
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end


