function sigma = Sqw_scatt_xs(s, Ei, M)
% Sqw_scatt_xs(s,Ei): compute the total scattering cross section for
%   incoming neutron energy. The S(|q|,w) should be the non-classical
%   dynamic structure factor. 
%
%   Such data sets are obtained from e.g. xray and neutron scattering 
%   experiments on isotropic density materials (liquids, powders, amorphous
%   systems). 
%
%   Data sets from analytical models and molecular dynamics simulations must 
%   be symmetrised in energy, and the detailed balance must be applied to 
%   take into account the material temperature on the inelastic part.
%
%   The S(q,w) is a dynamic structure factor aka scattering function.
%
%   The incident neutron energy is given in [meV], and may be computed:
%     Ei = 2.0721*Ki^2 = 81.8042/lambda^2 with Ki in [Angs-1] and lambda in [Angs]
%     
%   The S(|q|,w) is first restricted to the achievable dynamic range:
%     Ef         = Ei - w                                is positive
%     cos(theta) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf)  is within [-1:1]
%   and then integrated as XS = 1/2Ki^2 \int q S(q,w) dq dw
%
%   The computed value must then be multiplied by the tabulated bound cross-section
%   from e.g. Sears, Neut. News 3 (1992) 26.
%
%   When the weight M of the scattering unit is given, it is used to multiply the
%    cross section by the estimated Debye-Waller-like factor so that it equals
%    A/(A+1)]^2 at Ei=1eV to gradually go from the bound (thermal) to the free 
%    cross section at 1 eV. The threshold of 1eV is used by e.g. OpenMC.
%     W  = 2e-3*(log(M)-log(M+1));
%     DW = exp(W*Ei) = exp(2e-3*(log(M)-log(M+1))*Ei)
%   Above the ethermal energy threshold, the DW factor is kept fixed. 
%   For a poly-atomic scatterer, the effective mass is computed by weighting
%   with the bound cross sections for elements A, e.g.
%     r = sqrt(sum((A/(A+1))^2 * sigma_bound(A)) / sum(sigma_bound(A)));
%     M = r/(1-r)
%   For instance, for H2O (twice A=1 sigma=80.2 and A=18 sigma=4.2):
%     r = sqrt((2*(1/(1+1))^2*80.2+(18/(18+1))^2*4.2)/(2*80.2+4.2)) = 0.52
%     M = r/(1-r) = 1.06 i.e. scattering is mostly the hydrogen one.
%   WARNING: this factor is NOT the Debye-Waller factor exp(-<u2>Q2) !
%
% A classical S(|q|,w) obeys S(|q|,w) = S(|q|,-w) and is usually given
% on the positive energy side (w>0).
% The non classical S(q,w), needed by this function, can be obtained from a 
% classical S(q,w) (which is symmetric in energy) with e.g.:
%   extend to +/- energy range
%     s = Sqw_symmetrize(s); 
%   apply detailed balance (Bose factor). Omit T if you do not know it.
%     s = Sqw_Bosify(s, T);
%
% The positive energy values in the S(q,w) map correspond to Stokes processes, 
% i.e. material gains energy, and neutrons loose energy when scattered.
%
% input:
%   s: Sqw data set (non classical, with T Bose factor e.g from experiment)
%        e.g. 2D data set with w as 1st axis (rows, meV), q as 2nd axis (Angs-1).
%   Ei: incoming neutron energy [meV]
%   M: molar weight of the atom/molecule in [g/mol].
%     when given empty, it is searched 'weight' or 'mass' is the object.
%     Default is set to 0, i.e. the Debye-Waller factor is not taken into account.
% output:
%   sigma: cross section per scattering unit (scalar or iData)
%          to be multiplied afterwards by the bound cross section [barn]
%
% Example: sigma = Sqw_scatt_xs(s, 14.6)
%
% See also: Sqw_Bosify, Sqw_deBosify, Sqw_symmetrize, Sqw_dynamic_range
% (c) E.Farhi, ILL. License: EUPL.
  
  sigma = [];
  if nargin == 0, return; end
  if ~isa(s, 'iData'), s=iData(s); end

  % first look at the given arguments
  if nargin < 2,  Ei = []; end
  if nargin < 3,  M = []; end
  if isempty(Ei), Ei=14; end

  sigma = [];

  % handle array of objects
  if numel(s) > 1
    for index=1:numel(s)
      sigma = [ sigma feval(mfilename, s(index), Ei, M) ];
      s(index) = iData;
    end
    return
  end
  
  % get the S(q,w) on a meshgrid (and possibly rebin/interpolate), do checks
  s = meshgrid(s);
  s = Sqw_check(s);
  if isempty(s), return; end
  if isfield(s, 'classical') && s.classical == 1
    disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' seems to be "classical". You should apply s=Sqw_Bosify(s, temperature) before.' ]);
  end
  
  if isempty(Sqw_getT(s))
    disp([ mfilename ': WARNING: Temperature undefined: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' does not seem to have a Temperature defined. You may apply s=Sqw_Bosify(s, temperature) before.' ]);
  end
  if isempty(M) && isfield(s, 'weight')
    M       = s.weight;               % mass
  end
  
  w = s{1};
  
  if min(w(:)) * max(w(:)) >= 0
    disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' seems to only be defined on w<0 or w>0. You should apply  and then s=Sqw_Bosify(Sqw_symmetrize(s), temperature) before.' ]);
  end

  % do the total XS computation there...
  if numel(Ei) > 1
    sigma = [];
    if numel(Ei) == 2, Ei=logspace(log10(Ei(1)),log10(Ei(2)),20); end
    % loop for each incoming neutron energy
    for ie=1:numel(Ei)
      sigma = [ sigma Sqw_scatt_xs_single(s, Ei(ie), M) ];
    end
    sigma = iData(Ei, sigma);
    signa.Title = [ 'XS(' s.Title ')' ];
    title(sigma, 'XS [/barn/scatterer]');
    label(sigma, 'Signal', [ 'Total XS(' s.Title ')' ]);
    sigma.Label='XS';
    return
  else
    sigma = Sqw_scatt_xs_single(s, Ei, M);
  end
  
% ----------------------------------------------------------------------------
function sigma = Sqw_scatt_xs_single(s, Ei, M)

  % restrict to dynamic range
  s = Sqw_dynamic_range(s, Ei, [], 'checked');
  sq= trapz(s); % integrate over energy w {1} in the 2D Sqw

  % constants
  SE2V = 437.393377;        % Convert sqrt(E)[meV] to v[m/s]
  V2K  = 1.58825361e-3;     % Convert v[m/s] to k[1/AA]
  K2V  = 1/V2K;
  VS2E = 5.22703725e-6;     % Convert (v[m/s])**2 to E[meV]

  % compute int q.S(q)*2/Ki^2
  Ki    = SE2V*V2K*sqrt(Ei);
  q     = sq{1};
  sigma = trapz(q.*sq)/2/Ki^2; % integrate over q {1}
  
  % add the Debye-Waller factor computed to sweep bound to free cross section
  % at Ei=1 eV. Factor W is so that exp(-WE) = [A/(A+1)]^2 at threshold E=1 eV
  % DW factor is in q^2, i.e. in alpha i.e. in E, but this is not a Debye-Waller factor !
  if ~isempty(M) && M > 0
    Ei_threshold = 1000; % 1 eV
    W  = 2/Ei_threshold*(log(M)-log(M+1));
    Ei = min([ Ei Ei_threshold ]);
    DW = exp(W*Ei);
    sigma = sigma * DW;
  end
  
