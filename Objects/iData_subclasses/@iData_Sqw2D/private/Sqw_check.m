function s = Sqw_check(s, mode)
% Sqw_check: check if a 2D iData is a S(q,w).
%
% This routine can also convert automatically an input 
%         S(phi,t)  into S(phi,w)
%   and   S(phi,w)  into S(q,  w).
%
% conventions:
% omega = Ei-Ef = energy lost by the neutron
%    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%    omega < 0, neutron gains energy, anti-Stokes
%
% input:
%   s: Sqw data set
%        e.g. 2D data set with w as 1st axis (rows), q as 2nd axis.
%   mode: requested type of data. default is 'qw', but can also be 'ab' for S(alpha,beta)
%
% Example: sqw=Sqw_check('SQW_coh_lGe.nc');
% (c) E.Farhi, ILL. License: EUPL.

  if nargin == 0, return; end
  if ~isa(s, 'iData'), s = iData(s); end
  if nargin < 2, mode = []; end
  if isempty(mode), mode='qw'; end
  
  % handle array of objects
  if numel(s) > 1
    sqw = [];
    for index=1:numel(s)
      sqw = [ sqw feval(mfilename, s(index), mode) ];
    end
    s(index)=iData; % free memory
    s = sqw;
    return
  end
  
  if isempty(s), return; end
  
  % check if the data set is Sqw (2D)
  w_present=0;
  q_present=0;
  a_present=0;
  t_present=0;
  alpha_present=0;
  beta_present =0;
  if isa(s, 'iData') && ndims(s) >= 2
    for index=1:2
      lab = lower(label(s,index));
      def = getaxis(s, num2str(index));
      if ischar(def), lab = [ def ' ' lab ]; end
      if isempty(lab), lab=lower(getaxis(s, num2str(index))); end
      lab = strread(lab, '%s'); % split string into cell
      if strcmpm(lab, {'alpha','a'}) % strcmpm = multiple strcmpi is private below
        alpha_present=index;
      elseif strcmpm(lab, {'beta','b'})
        beta_present=index;
      elseif strcmpm(lab, {'wavevector','momentum','q','k','angs'})
        q_present=index;
      elseif strcmpm(lab, {'energy','frequency','w','e','mev','omega'})
        w_present=index;
      elseif strcmpm(lab, {'time','sec','t','tof','channels','channel'})
        t_present=index;
      elseif strcmpm(lab, {'angle','deg','theta','phi','2theta','2thetha'})
        a_present=index;
      end
    end
  end
  
  % search for Sqw parameters for further conversions
  if ~isfield(s, 'parameters')
    [s, parameters] = Sqw_parameters(s);
  end
  % conversions: S(a,b) -> S(q,w) when expecting 'qw'
  if alpha_present && beta_present && ~w_present && ~q_present
    if strcmp(mode, 'qw')
      disp([ mfilename ': S(alpha,beta) data set detected: converting to S(q,w)' ]);
      s = Sab2Sqw(s); % convert from S(alpha,beta) to S(q,w)
      s = feval(mfilename, s, mode);
      return
    end
  end
  % conversions: S(q,w) -> S(a,b) when expecting 'ab'
  if ~alpha_present && ~beta_present && w_present && q_present
    if strcmp(mode, 'ab')
      disp([ mfilename ': S(q,w) data set detected: converting to S(alpha,beta)' ]);
      s = Sqw2Sab(s); % convert from S(alpha,beta) to S(q,w)
      s = feval(mfilename, s, mode);
      return
    end
  end
  if ~w_present && t_present
    % convert from S(xx,t) to S(xx,w): t2e requires L2=Distance
    disp([ mfilename ': time/channel data set detected: converting axis "' label(s,t_present) '" to energy.' ]);
    s = Sqw_t2e(s, t_present);
    if ~isempty(s), w_present = t_present; end
  end
  if ~q_present && a_present && w_present
    % convert from S(phi,w) to S(q,w). WARNING: phi = 2*theta is the angle at the detector.
    disp([ mfilename ': S(phi,w) data set detected: converting angle axis "' label(s,a_present) '" to wavevector.' ]);
    s = Sqw_phi2q(s, [], a_present, w_present);
    if ~isempty(s), q_present = a_present; end
  end
  if (strcmp(mode, 'qw') && (~w_present || ~q_present)) ...
  || (strcmp(mode, 'ab') && (~alpha_present || ~beta_present))
    disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
    disp('    does not seem to be an isotropic neutron scattering law 2D object. Ignoring.');
    disp([ '    Energy    axis present:' num2str(w_present) ])
    disp([ '    Time      axis present:' num2str(t_present) ])
    disp([ '    Angle     axis present:' num2str(a_present) ])
    disp([ '    Momentum  axis present:' num2str(q_present) ])
    disp([ '    Alpha     axis present:' num2str(alpha_present) ])
    disp([ '    Beta      axis present:' num2str(beta_present) ])
    s = [];
    return
  end

  % check if we need to transpose the S(q,w)
  if (w_present==2 && q_present==1) || (beta_present==2 && alpha_present==1)
    s = transpose(s);
  end
  
  % this is the weighting for valid data
  S.type='()';
  S.subs={ ~isfinite(s) };
  s = subsasgn(s, S, 0);
  
  % check 'classical' and 'symmetric'
  if isfield(s,'classical') 
    classical0 = get(s,'classical');
  else
    classical0 = [];
  end
  
  % handle ENDF keywords
  if isfield(s,'LASYM')
    setalias(s,'classical', ~get(s,'LASYM'));
    classical0  = get(s,'classical');
  end
  if numel(classical0) > 1, classical0=classical0(1); end
  
  w  = getaxis(s,1);
  % checks that we have 0<w<Ei for Stokes, and w<0 can be lower (anti-Stokes)
  if any(w(:) < 0) && any(w(:) > 0) && strcmp(mode, 'qw')
    % for experimental data we should get w1=max(w) < Ei and w2 can be as low as
    % measured (neutron gains energy from sample).
    
    
    w1 = max(w(:)); w2 = max(-w(:)); % should have w1 < w2
    if w1 > w2*2
      % we assume the measurement range is at least [-2*Ei:Ei]
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
      disp('    indicates that the energy range is mostly in the positive side.')
      disp('    Check that it corresponds with the neutron loss/sample gain Stokes side');
      disp('    and if not, revert energy axis with e.g. setaxis(s, 1, -s{1})');
    end
  end

  % can we guess if this is classical data ? get temperature ?
  w = getaxis(s,1);

  if any(w(:) < 0) && any(w(:) > 0)  && strcmp(mode, 'qw')
    % compute Temperature from the detailed balance (Bose)
    T = Sqw_getT(s, [], 'Bose');
    if isempty(T), T=nan; end

    if isfinite(T) && T < -0.1
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
      disp([ '    indicates a negative temperature T=' num2str(T) ' K. ' ]);
      disp(  '    Check the definition of the energy: Stokes=neutron losses energy=positive energy side');
      T = Sqw_getT(s);
    end

    % temperature stored ?
    T0        = Sqw_getT(s);

    % log_s_ratio should be about 0 if S(q,w) is symmetric
    classical = [];
    if ~isfinite(T) && ~isnan(T)
      classical = 1;
      T         = Sqw_getT(s);
    elseif T <= 0.1
      classical = 0;
    elseif T > 3000
      classical = 1;
    else 
    end

    % display warnings when mismatch is found
    if ~isempty(classical0) && ~isempty(classical) && classical0 ~= classical
      if   classical0, classical_str='classical/symmetric';
      else             classical_str='experimental/Bose/quantum/asymmetric'; end
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
      disp(['    indicates a ' classical_str ' S(|q|,w) 2D object, but the analysis of the data shows it is not.' ]);
    elseif isempty(classical0) && ~isempty(classical)
      setalias(s,'classical', classical);
    end

    if ~isempty(T0) && ~isempty(T) && ~isnan(T) && T>0 && ~(0.9 < T/T0 & T/T0 < 1.1)
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' S(|q|,w) 2D object from ' s.Source ]);
      disp(['    indicates a Temperature T=' num2str(T0) ' [K], but the analysis of the data provides T=' num2str(T) ' [K].' ]);
    end
    if isempty(T0) && ~isempty(T) && ~isnan(T) && T > 0.1 && T < 3000
      disp([ mfilename ': INFO: Setting temperature T=' num2str(T) ' [K] for data set ' s.Tag ' ' s.Title ' S(|q|,w) 2D object from ' s.Source ]);
      s.Temperature = T;
    end
    
  end % energy axis has +/- 
  
  if isfield(s,'classical') 
    classical = get(s,'classical');
    if numel(classical) > 1, s = setalias(s, 'classical', classical(1)); end
  end

% ------------------------------------------------------------------------------
function flag=strcmpm(str, words)
% multiple strcmp
%
% input:
%   str:   string or cellstr of tokens
%   words: string or cellstr of words to search

  flag = false;
  if ischar(str),   str=strread(str, '%s','delimiter',' ,; $()[]{}=|<>&"/\:"'''); end
  if ischar(words), words = strread(words, '%s'); end;
  for index=1:numel(words)
    if any(strcmpi(str, words{index}))
      flag = true; return;
    end
  end
