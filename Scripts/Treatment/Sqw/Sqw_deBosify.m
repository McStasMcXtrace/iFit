function s = Sqw_deBosify(s, T)
% Sqw_deBosify: remove Bose factor (detailed balance) from an 'experimental' data set.
%  In principle the resulting data set is 'classical' that is S(q,w) = S(q,-w)
%
% input:
%   s: Sqw data set (non classical, with T Bose factor e.g from experiment)
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

  % check if the data set is Sqw (2D)
  w_present=0;
  q_present=0;
  if isa(s, 'iData') && ndims(s) == 2
    for index=1:2
      lab = lower(label(s,index));
      if any(strfind(lab, 'wavevector')) || any(strfind(lab, 'q')) || any(strfind(lab, 'Angs'))
        q_present=index;
      elseif any(strfind(lab, 'energy')) || any(strfind(lab, 'w')) || any(strfind(lab, 'meV'))
        w_present=index;
      end
    end
  end
  if ~w_present || ~q_present
    disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' does not seem to be an isotropic S(|q|,w) 2D object. Ignoring.' ]);
    return
  end

  if isempty(T),  T = Sqw_getT(s); end

  % test if classical
  if ~isempty(findfield(s, 'classical'))
    if s.classical == 1
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' seems to already be classical. The detailed balance removal may be wrong.' ]);
    end
  end

  % check if we need to transpose the S(q,w)
  if w_present==2 && q_present==1
    s = transpose(s);
  end
  
  % get symmetric from experimental data
  s = Sqw_Bosify(s, -T);
  s.Temperature = T;
  s.classical = 1;
  
