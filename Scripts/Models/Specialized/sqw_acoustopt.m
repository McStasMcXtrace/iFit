function signal=sqw_acoustopt(varargin)
% model = sqw_acoustopt(p, h,k,l,w, {signal}) : acoutsic/optic dispersion(HKL) with DHO(energy)
%
%   iFunc/sqw_acoustopt: a 4D S(q,w) with a 3D HKL dispersion with 
%      quadratic dependency, and a DHO line shape. 
%      This dispersion corresponds with a local description of an excitation, 
%      with its minimum around an (H0,K0,L0,E0) point.
%
%   The model requires to define a direction corresponding with a Slope1 dependency
%   as well as a second direction. An ortho-normal coordinate basis is then derived. 
%   All HKL coordinates are in rlu, and energies are in meV.
%
%   The dispersion has the form:
%      w(q) = sqrt(E0^2 + Slope^2*(q-HKL0)^2)
%   so that when the dispersion is linear for E0=0 or far from HKL0, else it is 
%   quadratic.
%
%   To define a pure acoustic mode, use (and adjust slopes):
%      sqw_acoustopt(0)                                 % minimum is E0=0 at q=0
%   To define an optical mode with energy E0 at Q=0, use:
%      sqw_acoustopt([ E0 ])                            % minimum is E0 at q=0
%   To define a mode which has its minimum E0 at a given HKL location, use:
%      sqw_acoustopt([ H K L E0 ])
%
%   When creating the Model, the following syntax is possible:
%      sqw_acoustopt(E0)           centers the excitation at q=0 with energy E0
%      sqw_acoustopt([ h k l E0 ]) centers the excitation at q=[H K L] and energy E0
%   You can of course tune other parameters once the model object has been created.
%
% WARNING: Single intensity and line width parameters are used here.
%
% To model more than one branch, just add these models together.
%
% Example:
%   s=sqw_acoustopt(5); qh=linspace(0,.5,50);qk=qh; ql=qh'; w=linspace(0.01,10,50);
%   f=iData(s,s.p,qh,qk,ql,w); plot3(log(f(:,1,:,:)));
%
% Reference: https://en.wikipedia.org/wiki/Phonon
%
% input:  p: sqw_acoustopt model parameters (double)
%             p(1) = DC_Hdir1         Slope1 dispersion direction, H [rlu]
%             p(2) = DC_Kdir1         Slope1 dispersion direction, K [rlu]
%             p(3) = DC_Ldir1         Slope1 dispersion direction, L [rlu]
%             p(4) = DC_Hdir2         Slope2 dispersion direction, H (transverse) [rlu]
%             p(5) = DC_Kdir2         Slope2 dispersion direction, K (transverse) [rlu]
%             p(6) = DC_Ldir2         Slope2 dispersion direction, L (transverse) [rlu]
%             p(7) = DC_Slope1        Dispersion slope along 1st axis [meV/rlu]
%             p(8) = DC_Slope2        Dispersion slope along 2nd axis (transverse) [meV/rlu]
%             p(9) = DC_Slope3        Dispersion slope along 3rd axis (vertical) [meV/rlu]
%             p(10)= Ex_H0            Minimum of the dispersion, H [rlu]
%             p(11)= Ex_K0            Minimum of the dispersion, K [rlu]
%             p(12)= Ex_L0            Minimum of the dispersion, L [rlu]
%             p(13)= Ex_E0_Center     Minimum of the dispersion, Energy [meV]
%             p(14)= DHO_Amplitude
%             p(15)= DHO_Damping      Excitation damping, half-width [meV]
%             p(16)= DHO_Temperature  Temperature [K]
%             p(17)= Background   
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

signal.Name           = [ 'Sqw_acoustopt acoustic/optic dispersion(HKL) with DHO line shape [' mfilename ']' ];
signal.Description    = 'A 4D S(q,w) with a 3D HKL quadratic dispersion, and a DHO line shape. This dispersion corresponds with a local description of an excitation, centered around an (H0,K0,L0,E0) point.';

signal.Parameters     = {  ...
'DC_Hdir1         Slope1 dispersion direction, H [rlu]' ...
'DC_Kdir1         Slope1 dispersion direction, K [rlu]' ...
'DC_Ldir1         Slope1 dispersion direction, L [rlu]' ...
'DC_Hdir2         Slope2 dispersion direction, H (transverse) [rlu]' ...
'DC_Kdir2         Slope2 dispersion direction, K (transverse) [rlu]' ...
'DC_Ldir2         Slope2 dispersion direction, L (transverse) [rlu]' ...
'DC_Slope1        Dispersion slope along 1st axis [meV/rlu]' ...
'DC_Slope2        Dispersion slope along 2nd axis (transverse to 1st, in plane) [meV/rlu]' ...
'DC_Slope3        Dispersion slope along 3rd axis (transverse to 1st, vertical) [meV/rlu]' ...
'Ex_H0            Minimum of the dispersion, H [rlu]' ...
'Ex_K0            Minimum of the dispersion, K [rlu]' ...
'Ex_L0            Minimum of the dispersion, L [rlu]' ...
'Ex_E0_Center     Minimum of the dispersion, Energy [meV]' ...
'DHO_Amplitude' ...
'DHO_Damping      Excitation damping, half-width [meV]' ...
'DHO_Temperature  Temperature [K]' ...
'Background' };
  
signal.Dimension      = 4;         % dimensionality of input space (axes) and result

signal.Expression     = { ...
'% define Bragg point, difference from Phonon pos', ...
'U1 = p(1:3);    % dir1 direction', ...
'U2 = p(4:6);    % dir2 direction', ...
'HKLEph = p(10:13);  % center of measurement (phonon)', ...
'% create an ortho-normal basis from lin and quad directions', ...
'U3=cross(U1,U2); % vertical', ...
'U2=cross(U3,U1);', ...
'U3=U3/norm(U3);', ...
'U2=U2/norm(U2);', ...
'U1=U1/norm(U1);', ...
'% now compute the dispersion terms b[xyz]=Bragg location, d[xyz]=excitation shift/Bragg', ...
'bx = round(x); dx = (x - bx - HKLEph(1));  % [rlu]', ...
'by = round(y); dy = (y - by - HKLEph(2));', ...
'bz = round(z); dz = (z - bz - HKLEph(3));', ...
'% dispersion terms: linear and quadratic', ...
'dvx = (dx*U1(1) +dy*U1(2) +dz*U1(3) ).^2*p(7)^2*sign(p(7)); % meV2, quad long', ...
'dvy = (dx*U2(1) +dy*U2(2) +dz*U2(3) ).^2*p(8)^2*sign(p(8)); % meV2, quad transv', ...
'dvz = (dx*U3(1) +dy*U3(2) +dz*U3(3) ).^2*p(9)^2*sign(p(9)); % meV2, quad vert', ...
'% handle curvature sign', ...
'wqs = HKLEph(4)^2 + dvx + dvy + dvz;', ...
'wq  = real(sqrt(wqs));', ...
'T   = p(16);', ...
'this.UserData.maxFreq=max(wq(:));', ...
'if T<=0, T=300; end', ...
'n_w    = 1./(exp(t*11.609/T)-1);', ...
'signal = p(17)+p(14)*p(15)*t.*(n_w+1)./((t.^2-wqs).^2+(p(15)*t).^2);' };

signal=iFunc(signal);

p = [ 1 0 0 ...
      0 1 0 ...
      20 10 10 ...
      0 0 0 5 ...
      1 0.05 10 0];
if nargin == 0
  signal.ParameterValues=p;
elseif nargin == 1 && isnumeric(varargin{1})
  if length(varargin{1}) == 1 % [ Emax ]
    p(13) = varargin{1}; 
  elseif length(varargin{1}) == 4 % [HKLE]
    p(10:13) = varargin{1};
    % [ 0 0 0 E0 ] optical
    % [ h k l E0 ] acoustic
  end
  signal.ParameterValues = p;
elseif nargin > 1
  signal = signal(varargin{:});
end

