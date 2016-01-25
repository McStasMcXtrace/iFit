function s = Sqw_deBosify(s, T)
% Sqw_deBosify: remove Bose factor (detailed balance) from an 'experimental' data set.
%  In principle the resulting data set is 'classical' that is S(q,w) = S(q,-w)
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
%  Bose factor: n(w) = 1./(exp(w*11.605/T) -1) ~ exp(-w*11.605/T)
%               w in [meV], T in [K]
%
% Example: s = Sqw_deBosify(s, 300);
%
% See also: Sqw_Bosify, Sqw_symmetrize, Sqw_dynamic_range, Sqw_total_xs

  if nargin == 1, T = []; end

  % handle array of objects
  if numel(s) > 1
    sqw = [];
    for index=1:numel(s)
      sqw = [ sqw feval(mfilename, s(index), T) ];
    end
    s = sqw;
    return
  end

  s = Sqw_check(s);
  if isempty(s), return; end

  % test if classical
  if isfield(s,'classical') || ~isempty(findfield(s, 'classical'))
    if s.classical == 1
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' seems to already be classical. The detailed balance removal may be wrong.' ]);
    end
  end
  
  % get symmetric from experimental data
  s = Sqw_Bosify(s, -T);
  setalias(s, 'Temperature', T);
  setalias(s, 'classical', 1);
  
