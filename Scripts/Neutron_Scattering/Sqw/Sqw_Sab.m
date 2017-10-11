function Sab = Sqw_Sab(s, M, T)
% Sqw = Sqw_Sab(Sqw, M, T)
%  Sqw_Sab: convert an S(alpha,beta) into a 2D S(q,w). 
%
%  The S(alpha,beta) is a representation of the dynamic structure factor 
%  using unitless momentum and energy variables defined as:
%     alpha= h2q2/2MkT  = (Ei+Ef-2*mu*sqrt(Ei*Ef))/AkT
%     beta = -hw/kT     = (Ef-Ei)/kT
%     A    = M/m
%     mu   = cos(theta) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf)
%  
% input:
%   s:  S(alpha,beta) data set e.g. 2D data set with alpha as 1st axis (rows), beta as 2nd axis.
%   M:  molar weight of the atom/molecule in [g/mol].
%     when omitted or empty, it is searched as 'weight' or 'mass' is the object.
%   T: when given, Temperature to use. When not given or empty, the Temperature
%      is searched in the object. The temperature is in [K]. 1 meV=11.605 K.
% output:
%   Sqw: S(q,w) 2D data set, classical/symmetric (iData)
%
% conventions:
% w = omega = Ei-Ef = energy lost by the neutron [meV]
%    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%    omega < 0, neutron gains energy, anti-Stokes
%
% references: M. Mattes and J. Keinert, IAEA INDC(NDS)-0470, 2005.
%             R. E. MacFarlane, LA-12639-MS (ENDF-356), 1994.
%
% Example: Sqw = iData( fullfile(ifitpath,'Data','SQW_coh_lGe.nc') );
%          Sab = Sqw_Sab(Sqw,72.6,300);
%          plot(log(Sab)
%
% (c) E.Farhi, ILL. License: EUPL.

  Sab = [];
  if nargin == 0, return; end
  if ~isa(s, 'iData'), s=iData(s); end
  
  if nargin < 2, M=[]; end
  if nargin < 3, T=[]; end
  
  % handle array of objects
  if numel(s) > 1
    for index=1:numel(s)
      Sab = [ Sab feval(mfilename, s(index), M, T) ];
    end
    return
  end
  
  s = Sqw_check(s);
  if isempty(s), return; end
  
  if isempty(M) && isfield(s, 'weight')
    M       = s.weight;               % mass
  end
  if isempty(T)
    T = Sqw_getT(s);
  end
  
  % check input parameters
  if isempty(M)
    disp([ mfilename ': ERROR: Mass undefined: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' does not provide information about the material molar weight. Use Sqw_Sab(Sqw, M).' ]);
    return
  end
  if isempty(T) || T<=0
    disp([ mfilename ': WARNING: Temperature undefined: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' does not have any temperature defined. Using T=300K. Use Sqw_Sab(Sqw, M, T) for other setting.' ]);
    T = 300;
  end
  
  % test if classical
  if isfield(s,'classical') || ~isempty(findfield(s, 'classical'))
    if s.classical == 0
      fprintf(1, '%s: %s: Converting to classical/symmetric for T=%g [K] (Sqw_deBosify)\n', mfilename, s.Title, T);
      s = Sqw_deBosify(s, T); % make it classical/symmetric
    end
  end
  
  E = s{1}; 
  q = s{2}; 
  Z = s{0};                 % S(q,w) value (axis rank 0)
  
  if numel(q) ~= size(Z,2) || numel(E) ~= size(Z,1)
    s = meshgrid(s, 'vector');
    E = s{1}; 
    q = s{2};
  end

  % some constants 
  mn      = 1.675E-027;      % neutron mass [kg]
  kb      = 1.381E-023;      % Boltzmann [J/K]
  e       = 1.602E-019;      % [C]
  HBAR    = 1.05457168e-34;  % Plank/2PI

  fprintf(1, '%s: %s: q=[%g:%g] w=[%g:%g]\n', mfilename, s.Title, min(q(:)), max(q(:)), min(E(:)), max(E(:)));
  
  % constant 
  q2toE   = HBAR*HBAR/2/mn/e*1000*1e20; % = 2.072 = [Angs^-2] to [meV] 
  % the definition of alpha states that it is per M(atom/molecule), so we divide
  q2toE   = q2toE/M;
  C       = e/1000/kb/T;

  % transform the axis
  % Q2 = ki2 + kf2 − 2ki kf cos (2 theta)
  % S(a,b) da db = S(q,w) dq dw
  % S(a,b)       = S(q,w) * J avec J=(dq.dw)/(dalpha.dbeta)=1/q*cte
  
  for i=1:size(Z,1)   % E
    % Jacobian for S(q,w) -> S(a,b) is J=(dq.dw)/(dalpha.dbeta) = 1/(2*q*C^2*q2toE)
    Z(i,:)     = Z(i,:)./q;
  end
  
  % and finally multiply by the constants
  alpha = q.*q*C*q2toE;
  beta  = -E *C;
  Z     = Z/(2*C*C*q2toE);

  % create new data set, and display it
  Sab=iData(alpha,beta,Z);  % Z(alpha,beta)
  Sab.Title = [ 'Sab(' s.Title ')' ];
  Sab.Temperature= T;
  Sab.weight     = M;
  Sab.classical  = 1;
  
  
  Sab.Label='Sab';
  title(Sab, 'S(\alpha,\beta)');
  ylabel(Sab,'alpha [h2q2/2MkT]');
  xlabel(Sab,'beta [-hw/kT]');
  Sab = transpose(Sab);
  % transfer available information compatible with ENDF MF7 MT4
  

