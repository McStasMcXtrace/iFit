function sqw = Sqw(s, M, T)
% sqw = Sqw(Sab, M, T)
%  iData_Sab: Sqw: convert a 2D S(alpha,beta) into an S(q,w).
%
%  The S(alpha,beta) is a representation of the dynamic structure factor 
%  using unitless momentum and energy variables defined as:
%     alpha= h2q2/2MkT  = (Ei+Ef-2*mu*sqrt(Ei*Ef))/AkT
%     beta = -hw/kT     = (Ef-Ei)/kT
%     A    = M/m
%     mu   = cos(theta) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf)
%  
% input:
%   s:  S(alpha,beta) data set e.g. 2D data set with beta as 1st axis (rows), alpha as 2nd axis (columns).
%   M:  molar weight of the atom/molecule in [g/mol].
%     when omitted or empty, it is searched as 'weight' or 'mass' is the object.
%   T: when given, Temperature to use. When not given or empty, the Temperature
%      is searched in the object. The temperature is in [K]. 1 meV=11.605 K.
% output:
%   sqw: S(q,w) 2D data set (iData_Sqw2D)
%
% conventions:
% w = omega = Ei-Ef = energy lost by the neutron [meV]
%    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%    omega < 0, neutron gains energy, anti-Stokes
%
% references: M. Mattes and J. Keinert, IAEA INDC(NDS)-0470, 2005.
%             R. E. MacFarlane, LA-12639-MS (ENDF-356), 1994.
%
% Example: sab = iData_Sab('SQW_coh_lGe.nc');
%          sqw = Sqw(sab,72.6,1235);
%          sab2= Sab(sqw);
%          subplot(log([Sqw Sab Sqw2))
%
% (c) E.Farhi, ILL. License: EUPL.

  sqw = [];
  
  if nargin < 2, M=[]; end
  if nargin < 3, T=[]; end
  
  % handle array of objects
  if numel(s) > 1
    for index=1:numel(s)
      sqw = [ sqw feval(mfilename, s(index), M, T) ];
    end
    return
  end

  if isempty(s), return; end
  
  if isempty(M) || M<=0
    M = Sqw_getT(s, {'Masses','Molar_mass','Mass','Weight'});
  end
  if isempty(T)
    T = Sqw_getT(s);
  end
  
  % check input parameters
  if isempty(M)
    disp([ mfilename ': ERROR: Mass undefined: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' does not provide information about the material molar weight. Use Sqw(Sab, M).' ]);
    return
  end
  if isempty(T) || T<=0
    disp([ mfilename ': ERROR: Temperature undefined: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' does not have any temperature defined. Use Sqw(Sab, M, T).' ]);
    return
  end
  
  beta       = getaxis(s,1); 
  alpha      = getaxis(s,2); 
  Z=getaxis(s,0);                 % S(a,b) value (axis rank 0)
  
  if numel(alpha) ~= size(Z,2) || numel(beta) ~= size(Z,1)
    s = meshgrid(s, 'vector');
    beta       = getaxis(s,1); 
    alpha      = getaxis(s,2); 
  end
  
  % we must have s.classical = 1
  if (isfield(s,'classical') || ~isempty(findfield(s, 'classical'))) && get(s,'classical') == 0
    fprintf(1, '%s: %s: Data set is quantum/experimental. Use deBosify to obtain the classical/symmetric one.\n', mfilename, s.Title);
  end

  % some constants 
  mn      = 1.675E-027;      % neutron mass [kg]
  kb      = 1.381E-023;      % Boltzmann [J/K]
  e       = 1.602E-019;      % [C]
  HBAR    = 1.05457168e-34;  % Plank/2PI
  
  % constant 
  q2toE   = HBAR*HBAR/2/mn/e*1000*1e20; % = 2.072 = [Angs^-2] to [meV] 
  % the definition of alpha states that it is per M(atom/molecule), so we divide
  q2toE   = q2toE/M;
  C       = e/1000/kb/T;

  % transform the axis
  % Q2 = ki2 + kf2 − 2ki kf cos (2 theta)
  % S(a,b) da db = S(q,w) dq dw
  % S(a,b)       = S(q,w) * J avec J=(dq.dw)/(dalpha.dbeta)=1/q*cte
  
  q = sqrt(alpha/(C*q2toE));
  E = -beta/C;
  
  % fprintf(1, '%s: %s: q=[%g:%g] w=[%g:%g]\n', mfilename, s.Title, min(q(:)), max(q(:)), min(E(:)), max(E(:)));
  
  for i=1:size(Z,1)   % E
    % Jacobian for S(q,w) -> S(a,b) is J=(dq.dw)/(dalpha.dbeta) = 1/(2*q*C^2*q2toE)
    Z(i,:)     = Z(i,:).*q;
  end
  
  % and finally multiply by the constants
  
  Z     = Z*(2*C*C*q2toE);

  % create new data set, and display it
  sqw=iData(q,E,Z);  % Z(alpha,beta)
  sqw.Title = s.Title;
  setalias(sqw,'Temperature', T, '[K] Temperature');
  setalias(sqw,'weight',      M, '[g/mol] Material molar weight');
  if ~isfield(s,'classical')
    setalias(sqw,'classical',   1, '[0=from measurement, with Bose factor included, 1=from MD, symmetric]');
  else
    setalias(sqw,'classical',  get(s,'classical'), '[0=from measurement, with Bose factor included, 1=from MD, symmetric]');
  end
  sqw.Label='Sqw';
  title(sqw, 'S(q,w)');
  ylabel(sqw,'wavevector [Angs-1]');
  xlabel(sqw,'energy [meV]');
  sqw = transpose(sqw);
  
  % copy initial aliases and UserData
  sqw.UserData = s.UserData;
  f = getalias(s);
  for index=1:numel(getalias(s))
    if ~isfield(sqw, f{index})
      [link, lab] = getalias(s, f{index});
      sqw = setalias(sqw, f{index}, link, lab);
    end
  end
  % transfer available information compatible with ENDF
  sqw = iData_Sqw2D(sqw);
  
  if nargout == 0
    fig=figure; 
    h  =plot(sqw); 
    set(fig, 'NextPlot','new');
  end

