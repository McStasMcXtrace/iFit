function signal=sqw_sine3d(varargin)
% model = sqw_sine3d(p, h,k,l,w, {signal}) : sine dispersion(HKL) with DHO(energy)
%
%   iFunc/sqw_sine3d: a 4D S(q,w) with a 3D HKL sine dispersion, and a DHO line
%      shape. The sine dispersion can be tuned with energy gaps at zone centre
%      and boundary, as well as a periodicity of the sine wave.
%
% WARNING: Single intensity and line width parameters are used here.
%
% Along a principal Q axis, the dispersion has the form:
%       w = E0 + (E1-E0)*sin(Q_freq*pi*(Q-Q0));
% This is a sine wave which goes from w=E0 at Q=Q0, up to w=E1 at 1/2Q_freq.
% The inter-atomic distance between the scattering units (atoms) is thus:
%       d=a*Q_freq/2 [in Angs, with a=lattice parameter]
%
% A spin wave could for instance mostly use Q0=0, Q_freq=1,  E0=0, E1>0
%   For a simple ferromagnetic magnon the gap width is E1-E0=4JS 
%     with J=exchange energy and S=magnetic moment of spins.
% An acoustic branch could use           Q0=0, Q_freq=.5, E0=0, E1>0
%   The sound velocity in acoustic branches is:
%       c=E1*Q_freq*pi*1.519e2*a/2 [in m/s, with a=lattice parameter]
% An optical branch could use            Q0=0, Q_freq=.2, E0>E1 E1>0
%
% To shift the minimum/maximum Q of the dispersion, move QH0,QK0,QL0
% To change the extent of the dispersion in Q, vary QH_freq,QK_freq,QL_freq
% To change to minimum and maximum energy, move E0 and E1_qh,E1_qk,E1_ql
% To model an incomensurate dispersion, move both QH0,QK0,QL0 and the
%   QH_freq,QL_freq,QK_freq to incomensurate (non rational) values.
%
% To quickly create predefined models, use:
%   model=sqw_sine3d(Emax)      creates an acoustic dispersion up to Emax
%   model=sqw_sine3d([ E0 E1 ]) creates an optical dispersion from E0 to E1
%   model=sqw_sine3d([ E0 E1 Q_freq ])    creates a dispersion from E0 to E1
%                                      with given Q frequency, e.g. .5, 1 or 2
% Of course, as for any iFunc model, parameter values can be changed afterwards.
%
% To model more than one branch, just add these models together, e.g.:
%     disp3 = sqw_sine3d(5) + sqw_sine3d([ 10 8 ]) + sqw_sine3d([ 2 4 1 ])
% which is an acoustic plus optical branch, and a spin wave with 2 meV gap.
%
% Example:
%   s=sqw_sine3d(5); qh=linspace(0,1,50);qk=qh; ql=qh'; w=linspace(0.01,10,50);
%   f=iData(s,s.p,qh,qk,ql,w); plot(log(f(:,:,1,:)));
%
% Reference: https://en.wikipedia.org/wiki/Phonon
%
% input:  p: sqw_sine3d model parameters (double)
%             p(1)= E1_qh   energy at QH half period [meV]
%             p(2)= E1_qk   energy at QK half period [meV]
%             p(3)= E1_ql   energy at QL half period [meV]
%             p(4)= E0      zone-centre energy gap [meV]
%             p(5)= QH0     QH zone-centre [rlu]
%             p(6)= QK0     QK zone-centre [rlu]
%             p(7)= QL0     QL zone-centre [rlu]
%             p(8)= QH_freq QH frequency [multiples of pi]
%             p(9)= QK_freq QK frequency [multiples of pi]
%             p(10)=QL_freq QL frequency [multiples of pi]
%             p(11)=Gamma   dispersion DHO half-width in energy [meV]
%             p(12)=Temperature of the material [K]
%             p(13)=Amplitude
%             p(14)=Background (constant)
%          or p='guess'
%         qh: axis along QH in rlu (row,double)
%         qk: axis along QK in rlu (column,double)
%         ql: axis along QL in rlu (page,double)
%         w:  axis along energy in meV (double)
%    signal: when values are given, a guess of the parameters is performed (double)
% output: signal: model value
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot, gauss, sqw_phonons, sqw_cubic_monoatomic, sqw_vaks
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>
% (c) E.Farhi, ILL. License: EUPL.

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
signal = iFunc_Sqw4D(signal); % overload Sqw4D flavour

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

