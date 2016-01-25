function g = Sqw_gDOS(s,T)
% Sqw_gDOS: compute the generalised density of states (gDOS)
%  The gDOS is an approximation of the vibrational spectra (DOS).
%  This routine should be applied on an incoherent dynamic S(q,w) data set.
%
%  The S(q,w) is a dynamic structure factor aka scattering function.
%
% input:
%   s: Sqw data set (non classical, including T Bose factor e.g from experiment)
%        e.g. 2D data set with w as 1st axis (rows), q as 2nd axis.
%   T: when given, Temperature to use for Bose. When not given, the Temperature
%      is searched in the object.
%
% conventions:
% omega = Ei-Ef = energy lost by the neutron
%    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%    omega < 0, neutron gains energy, anti-Stokes
%
% references: Price J. et al, Non Cryst Sol 92 (1987) 153
%             Bellisent-Funel et al, J. Mol. Struct. 250 (1991) 213
%
% Example: s = Sqw_gDOS(s);


  if nargin == 1, T = []; end
  g = [];

  % handle array of objects
  if numel(s) > 1
    g = [];
    for index=1:numel(s)
      g = [ g feval(mfilename, s(index), T) ];
    end
    return
  end
  
  s = Sqw_check(s);
  if isempty(s), return; end

  if isempty(T),  T = Sqw_getT(s); end

  % test if classical
  if isfield(s,'classical') || ~isempty(findfield(s, 'classical'))
    if s.classical == 1
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' seems to be classical. The gDOS computation may be wrong. Apply Sqw_Bosify first.' ]);
    end
  end
  
  % check if we need to transpose the S(q,w)
  if w_present==2 && q_present==1
    s = transpose(s);
  end
  
  w = s{1};
  q = s{2};
  kT= T/11.6045;

  g = s.*(w.*(1 - exp(-w/kT)))./q.^2;  % from Price 1987
  
  % reset axes
  g{1}=w; g{2}=q;
  xlabel(g, xlabel(s));
  ylabel(g, ylabel(s));
  title(g,'gDOS');
  if isempty(g.Label), g.Label='gDOS'; end
  
  % this is the weighting for valid data
  g(g <= 0) = 0; g(~isfinite(g)) = 0;

  % compute integral
  %g = trapz(g,2); % on energy
  %g = g./trapz(g);
  setalias(g, 'Temperature', T);
  

