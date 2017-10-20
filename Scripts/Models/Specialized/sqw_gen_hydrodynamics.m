function signal=sqw_gen_hydrodynamics(varargin)
% model = sqw_gen_hydrodynamics(p, h,k,l,w, {signal}) : Generalized Hydrodynamics model
%
%  This is a 1D model to fit the parameters at a given 'q' value.
%
% Reference: U. Bafile, E. Guarini, and F. Barocchi, Phys. Rev. E 73, 061203
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot, gauss, sqw_phonons, sqw_cubic_monoatomic, sqw_vaks
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>
% (c) E.Farhi, ILL. License: EUPL.

% input parameters
sq;     % vector(q) structure factor
a0;     % vector(q) amplitude  of central Lorentzian, thermal non propagating
z0;     % vector(q) half-width of central Lorentzian, thermal non propagating
zs;     % vector(q) damping    of the sound mode
ws      % vector(q) frequency  of the sound mode
as;     % amplitude of the inelastic lines
bs;     % asymmetry factor of the inelastic lines
T;      % Temperature in [K]

% first search iq = index in q for the given 'x' momentum values
q = x;
w = y;

% Eq (4) from PRB 88(10):104201 · September 2013 
%             DOI: 10.1103/PhysRevB.88.104201
SGH(q,w)  = sq(iq)/pi.*( ...
              a0(iq) .* z0(iq) / (z0(iq).^2 + w.^2) ...
            + as .* (zs(iq)+bs*(w + ws(iq))) ./ (zs.^2 + (w + ws).^2) ...
            + as .* (zs(iq)-bs*(w - ws(iq))) ./ (zs.^2 + (w - ws).^2) );
         
beta      = w*11.608/T;  % hw/kT
bose      = beta*w./(1-exp(-w*beta)); % detailed balance (Bose factor)

S(q,w)  = bose .* SGH(q,w);




signal.Name           = [ 'Sqw_sine3d sine dispersion(HKL) with DHO line shape [' mfilename ']' ];
signal.Description    = 'A HKL sine wave dispersion with tunable energy gap at zone centre and boundary, and DHO line shape';

signal.Parameters     = {  ...
  'E1_qh energy at QH half period [meV]' ...
  'E1_qk energy at QK half period [meV]' ...
  'E1_ql energy at QL half period [meV]' ...
  'E0 zone-centre energy gap [meV]' ...
  'QH0 QH zone-centre [rlu]' ...
  'QK0 QK zone-centre [rlu]' ...
  'QL0 QL zone-centre [rlu]' ...
  'QH_freq QH frequency [multiples of pi]' ...
  'QK_freq QK frequency [multiples of pi]' ...
  'QL_freq QL frequency [multiples of pi]' ...
  'Gamma Damped Harmonic Oscillator width in energy [meV]' ...
  'Temperature [K]' ...
  'Amplitude' 'Background' };
  
signal.Dimension      = 4;         % dimensionality of input space (axes) and result

m1 = @(x,s) sum(s(:).*x(:))/sum(s(:));
m2 = @(x,s) sqrt(abs( sum(x(:).*x(:).*s(:))/sum(s(:)) - m1(x,s).^2 ));

signal.Guess          = @(x,y,z,t,signal)[ ...
  m2(t,signal) m2(t,signal) m2(t,signal) ...
  std(t(:)) round(mean(x(:))) round(mean(y(:))) round(mean(z(:))) ...
  1 1 1 std(t(:))/100 50 1 0 ];        % default parameters
  
% get code to read xyzt and build HKL list and convolve DHO line shapes
  [script_hkl, script_dho] = sqw_phonons_templates;

signal.Expression     = { ...
'% x=qh; y=qk; z=ql; t=w', ...
script_hkl{:}, ...
'wqx = (p(1)^2-p(4)^2)*sin(p(8) *pi*(x-p(5))).^2;', ...
'wqy = (p(2)^2-p(4)^2)*sin(p(9) *pi*(y-p(6))).^2;', ...
'wqz = (p(3)^2-p(4)^2)*sin(p(10)*pi*(z-p(7))).^2;', ...
'wq  = p(4)^2+wqx+wqy+wqz;', ...
'FREQ=wq(:);', ...
'clear wqx wqy wqz wq', ...
'Gamma=p(11); T=p(12); Amplitude=p(13); Bkg=p(14);', ...
script_dho{:} };

signal=iFunc(signal);

if nargin == 0
  signal.ParameterValues=[ 4 6 8 0  0 0 0  .5 .5 .5  0.04  50  1 0 ];
elseif nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1 % [ Emax ]
    p = [varargin{1} varargin{1} varargin{1} 0  0 0 0  .5 .5 .5  0.04  50  1 0];
  elseif length(varargin{1}) == 2 % [ Emin Emax ]
    p = varargin{1};
    p = [p(2) p(2) p(2) p(1)  0 0 0  ...
         .5 .5 .5  0.04  50  1 0];
  elseif length(varargin{1}) == 3 % [ Emin Emax Q_freq ]
    p = varargin{1};
    p = [p(2) p(2) p(2) p(1)  0 0 0  ...
         p(3) p(3) p(3)  0.04  50  1 0];
  else
    p = varargin{1};
  end
  signal.ParameterValues = p;
elseif nargin > 1
  signal = signal(varargin{:});
end

