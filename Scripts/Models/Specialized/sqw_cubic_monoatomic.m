function signal=sqw_cubic_monoatomic(varargin)
% model = sqw_cubic_monoatomic(p, h,k,l,w, {signal}) : cubic monoatomic dispersion(HKL) with DHO(energy)
%
%   iFunc/sqw_cubic_monoatomic: a 4D S(q,w) with a HKL dispersion, and a DHO line
%      shape, for a simple monoatomic system on a cubic lattice. Only the 3 
%      acoustic modes are computed.
%
% WARNING:
%      Single intensity and line width parameters are used here.
%      The HKL position is relative to the closest Bragg peak. 
%
% Example:
% s=sqw_cubic_monoatomic([2 3]); qh=linspace(0,.5,50);qk=qh; ql=qh; w=linspace(0.01,10,51);
% f=iData(s,[],qh,qk,ql,w); scatter3(log(f(:,:,1,:)),'filled');
%
% References: https://en.wikipedia.org/wiki/Phonon
%  E. Meisterhofer <http://lampx.tugraz.at/~hadley/ss1/phonons/scdos/sc.pdf> 2013
%  Kittel, Solid State Physics, Wiley, New York, 7th ed. (1996), pp 99-103.
%
% input:  p: sqw_cubic_monoatomic model parameters (double)
%           p(1)=C_ratio C1/C2 force constant ratio first/second neighbours.
%           p(2)=E0      sqrt(C1/m) energy [meV]
%           p(3)=Gamma   dispersion DHO half-width in energy [meV]
%           p(4)=Temperature of the material [K]
%           p(5)=Amplitude
%           p(6)=Background (constant)
%         qh: axis along QH in rlu (row,double)
%         qk: axis along QK in rlu (column,double)
%         ql: axis along QL in rlu (page,double)
%         w:  axis along energy in meV (double)
%    signal: when values are given, a guess of the parameters is performed (double)
% output: signal: model value
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot, gauss, sqw_phonons, sqw_sine3d, sqw_vaks
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>
% (c) E.Farhi, ILL. License: EUPL.

signal.Name           = [ 'Sqw_cubic dispersion(HKL) for cubic monoatomic crystal [' mfilename ']' ];
signal.Description    = 'A phonon dispersion(HKL) for a cubic monoatomic crystal with DHO(energy) shape. From the dynamical matrix (2 neighbors).';

signal.Parameters     = {  ...
  'C_ratio C1/C2 force constant ratio first/second neighbours' ...
  'E0      sqrt(C1/m) energy [meV]' ...
  'Gamma   Dampled Harmonic Oscillator width in energy [meV]' ...
  'Temperature [K]' ...
  'Amplitude' 'Background' };
  
signal.Dimension      = 4;           % dimensionality of input space (axes) and result

signal.Guess          = [1 1 .1 10 1 0];        % default parameters

% get code to read xyzt and build HKL list and convolve DHO line shapes
  [script_hkl, script_dho] = sqw_phonons_templates;

signal.Expression     = { ...
  '% get parameter values',...
  'c_quot=p(1); E0=p(2);', ...
  script_hkl{:}, ...
  'M = nan(3);', ...
  'FREQ=zeros(size(HKL,1),3);',...
  'for index=1:size(HKL,1); kx=HKL(index,1); ky=HKL(index,2); kz=HKL(index,3);', ...
  '  M(1, 1) = 2*(1 - cos(2*pi*kx)) +2/c_quot*(2- cos(2*pi*kx)*cos(2*pi*ky)-cos(2*pi*kx)*cos(2*pi*kz)) ;', ...
  '  M(1, 2) = 2/c_quot* sin(2*pi*kx)* sin(2*pi*ky) ;', ...
  '  M(1, 3) = 2/c_quot* sin(2*pi*kx)* sin(2*pi*kz) ;', ...
  '  M(2, 1) = M(1, 2) ;', ...
  '  M(2, 2) = 2*(1 - cos(2*pi*ky)) +2/c_quot*(2- cos(2*pi*kx)*cos(2*pi*ky)-cos(2*pi*ky)*cos(2*pi*kz)) ;', ...
  '  M(2, 3) = 2/c_quot* sin(2*pi*ky)* sin(2*pi*kz) ;', ...
  '  M(3, 1) = M(1, 3) ;', ...
  '  M(3, 2) = M(2, 3) ;', ...
  '  M(3, 3) = 2*(1 - cos(2*pi*kz)) +2/c_quot*(2-cos(2*pi*kx)*cos(2*pi*kz)-cos(2*pi*ky)*cos(2*pi*kz)) ;', ...
  '  [eigvectors,eigvalues ] = eig(M);',...
  '  eigvalues = sqrt(diag(eigvalues))*E0;',...
  '  [dummy,sorted]=max(abs(eigvectors''));',...
  '  eigvalues = eigvalues(sorted);',...
  '  FREQ(index,:) = eigvalues;', ...
  'end % index in kx,ky,kz', ...
  'Gamma=p(3); T=p(4); Amplitude=p(5); Bkg=p(6);', ...
  script_dho{:}, ...
};

signal.UserData.properties.spacegroup        = 'Cubic (230)';
signal.UserData.properties.spacegroup_number = 230;

signal=iFunc(signal);
signal = iFunc_Sqw4D(signal); % overload Sqw4D flavour

if nargin == 0
  varargin{1} = [2 3];
elseif numel(varargin) == 1 && isnumeric(varargin{1})
  p = varargin{1};
  if numel(p) == 2, p =[ p(:)' .1 10 1 0 ]; end 
  signal.ParameterValues = p;
elseif nargin > 1
  signal = signal(varargin{:});
end

