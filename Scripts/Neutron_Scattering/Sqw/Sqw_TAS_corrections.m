function [data,V, eff] = Sqw_TAS_corrections(data)
  % extract the q,w axes and compute kf, lambda_f, theta_A for corrections
  % (c) E.Farhi, ILL. License: EUPL.
  
  if nargin == 0, return; end
  if ~isa(s, 'iData'), s=iData(s); end
  
  q       = data{2};
  w       = data{1};
  
  % constants
  KI      = 3.84;
  DA      = 3.355;
  DM      = 3.355;
  Phi     = 2.54e-2;
  P       = 5;
  RA      = iData('/usr/local/lib/mcstas-2.0/data/HOPG.rfl');
  
  lambda_i= 2*pi/KI;
  EI      = 81.8040/lambda_i/lambda_i;
  
  EF      = EI-w;  %w = EF-EI
  lambda_f= sqrt(81.8040./EF);
  KF      = 2*pi./lambda_f;
  
  A1      = asin(pi/DM/KI);
  A5      = asin(pi/DA./KF);
  
  
  % resolution volume:   RA(KF)*KF.^3./tan(A5)
  ra      = double(interp(RA,{KF}));
  ra(find(KF > 3.698)) = 0.59885932;
  V       = ra.*KF.^3./tan(A5);
  
  % detector efficiency: 1-exp(-0.0591*P*lambda_f*Phi)
  eff     = 1-exp(-0.0591*P*lambda_f*Phi);
  % [ KF V eff ]
  
  % apply corrections
  data = data./V./eff;
