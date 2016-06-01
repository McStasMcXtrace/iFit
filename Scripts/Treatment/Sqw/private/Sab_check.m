function s = Sab_check(s)
% Sab_check: check if a 2D iData is a S(alpha,beta).
%
% conventions:
% w = omega = Ei-Ef = energy lost by the neutron
%    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%    omega < 0, neutron gains energy, anti-Stokes
%
% input:
%   s: Sqw data set
%        e.g. 2D data set with w as 1st axis (rows), q as 2nd axis.

  if isempty(s), return; end
  
  % check if the data set is Sqw (2D)
  alpha_present=0;
  beta_present=0;
  if isa(s, 'iData') && ndims(s) == 2
    for index=1:2
      lab = lower(label(s,index));
      if isempty(lab), lab=lower(getaxis(s, num2str(index))); end
      if any(strfind(lab, 'alpha')) || strcmp(strtok(lab), 'a')
        beta_present=index;
      elseif any(strfind(lab, 'beta')) || strcmp(strtok(lab), 'b')
        alpha_present=index;
      end
    end
  end
  if ~alpha_present || ~beta_present
    disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
    disp('    does not seem to be an isotropic S(alpha,beta) 2D object. Ignoring.');
    s = [];
    return
  end

  % check if we need to transpose the S(q,w)
  if alpha_present==2 && beta_present==1
    s = transpose(s);
  end
  
  % this is the weighting for valid data
  s(~isfinite(s)) = 0;
  
  % check 'classical' and 'symmetric'
  for f={'classical', 'symmetric'}
    if isfield(s,f{1}) || ~isempty(findfield(s, f{1}))
      if isfield(s,f{1}) classical0 = s.(f{1});
      else 
        classical0 = findfield(s, f{1});
        if iscell(classical0), classical0=classical0{1}; end
        classical0 = get(s, classical0);
      end
      if numel(classical0) > 1 || (~isempty(classical0) && ~isfield(s,'classical'))
        classical0  = classical0(1);
        s.classical = classical0;
        break
      end
    end
  end

% ------------------------------------------------------------------------------
