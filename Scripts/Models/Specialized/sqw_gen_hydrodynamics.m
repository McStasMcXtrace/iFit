function signal=sqw_gen_hydrodynamics(varargin)
% model = sqw_gen_hydrodynamics(p, w, {signal}) : Generalized Hydrodynamics model
%
%  iFunc/sqw_gen_hydrodynamics Generalized Hydrodynamics model with Rayleigh-Brillouin triplet
%    This model can be used to fit e.g. a spectrum with a central non dispersive 
%    quasi elastic line  (lorentizian shape) and two dispersive (e.g. acoustic)
%    Stokes and anti-Stokes lines. This model is suited for e.g dense fluids at 
%    small Q momentum values.
%
%  This is a 1D model to fit the parameters at a given 'q' momentum value.
%  The second frequency moment is obtained from the model parameters as:
%    <w^2> = zs^2 + ws^2 - a0 [ (z0-zs)^2 + ws^2 ] = kT q^2/[m S(q)]
%
%  The S(q,w) at one momentum value is computed in its symmetrized expression 
%  (so-called classical). To get the 'true' quantum S(q,w) definition, use e.g.
%    sqw = sqw_gen_hydrodynamics .* bose;
%  where the Temperature is then given in [x units]. If 'x' is an energy in [meV]
%  then the Temperature parameter is T[K]/11.6045
%
%  To add a 'background' use e.g.
%    sqw = sqw_gen_hydrodynamics + constant('Background');
% or
%    sqw = sqw_gen_hydrodynamics + strline; 
%
%  To convolute with a pseudo-Voigt or Gaussian instrumental line width, use e.g.
%    sqw = convn(sqw_gen_hydrodynamics, gauss)
% or
%    sqw = convn(sqw_gen_hydrodynamics, pseudovoigt)
%
% Reference: 
%  U. Bafile, E. Guarini, and F. Barocchi, Phys. Rev. E 73, 061203 (2006)
%  E. Guarini et al,                       Phys. Rev. B 88, 104201 (2013)
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot, gauss, sqw_phonons, sqw_cubic_monoatomic, sqw_vaks
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>
% (c) E.Farhi, ILL. License: EUPL.

signal.Name           = [ 'Sqw_gen_hydrodynamics Generalized Hydrodynamics - Rayleight-Brillouin triplet [' mfilename ']' ];
signal.Description    = 'Generalized Hydrodynamics model with a central line, and two dispersive acoutic lines, suited for e.g. dense fluids';

signal.Parameters     = {  ...
  'sq     structure factor S(q)' ...
  'a0     amplitude(q)  of central Lorentzian, thermal non propagating' ...
  'z0     half-width(q) of central Lorentzian, thermal non propagating' ...
  'zs     damping(q)    of the sound mode' ...
  'ws     frequency(q)  of the sound mode' ...
  'as     amplitude(q) of the inelastic line' ...
  'bs     asymmetry factor of the inelastic line' ...
};
  
signal.Dimension      = 1;         % dimensionality of input space (axes) and result

% we must provide a p(xx) expression to tell iFunc to validate xx parameters
% then we use 'deal' on the cell array to distribute all parameters
signal.Expression = { ...
  '% the energy axis' ...
  'w = x; bs = p(7);' ...
  'pc =num2cell(p); [sq a0 z0 zs ws as bs] = deal(pc{:});' ...
  'signal  = sq/pi.*( ...' ...
  '              a0 * z0 ./ (z0.^2 + w.^2) ...' ...
  '            + as * (zs+bs*(w + ws)) ./ (zs.^2 + (w + ws).^2) ...' ...
  '            + as * (zs-bs*(w - ws)) ./ (zs.^2 + (w - ws).^2) );' ...   
};

m1 = @(x,s) sum(s(:).*x(:))/sum(s(:));
m2 = @(x,s) sqrt(abs( sum(x(:).*x(:).*s(:))/sum(s(:)) - m1(x,s).^2 ));

signal.Guess          = @(x,signal)[ ...
  max(signal(:)) ...
  1 m2(x,signal)/10 ...
  m2(x,signal)/10 m2(x,signal) median(signal(:)) 1 30 ];        % default parameters

signal=iFunc(signal);

if nargin == 1 && isnumeric(varargin{1})
  y.ParameterValues = varargin{1};
elseif nargin > 1
  y = y(varargin{:});
end

