function M=moments(self, select)
% iFunc_Sqw2D: m=moments(sqw): compute Sqw moments/sum rules (harmonic frequencies)
%
%   m = moments(sqw)
%
%  Compute the structure factor (moment 0), recoil energy (moment 1) and the
%    collective, harmonic and mean energy transfer dispersions.
%
%  The result is given as an iFunc array with models sets (see below) vs momentum.
%
% Reference: 
%   Helmut Schober, Journal of Neutron Research 17 (2014) pp. 109
%   Lovesey, Theory of Neutron Scattering from Condensed Matter, Vol 1, p180 eq. 5.38 (w0)
%   J-P.Hansen and I.R.McDonald, Theory of simple liquids Academic Press New York 2006.
%
% syntax:
%   m = moments(sqw)
%   m = moments(sqw, type)
%
% input:
%   sqw:  Sqw 2D model with q as 1st axis (Angs-1), w as 2nd axis (meV).
%   type: to select which moment to return. Can be left empty for all, or a
%         cellstr with any of: (e.g. {'sq','wc'})
%     'sq'     S(q) = \int S(q,w) dw = <S(q,w)>                 structure factor [moment 0]
%     'recoil' Er   = \int w*S(q,w) dw = <wS(q,w)> = h2q2/2M       recoil energy [moment 1]
%     'wc'     Wc   = sqrt(2kT*Er/S(q))                    collective/isothermal dispersion
%     'wl'     Wl                                          harmonic/longitudinal excitation
%     'wq'     Wq   = 2q*sqrt(kT/S(q)/M)                               mean energy transfer
%     'm2'     M2   = <w2S(q,w)>                                                 [moment 2]
%     'm3'     M3   = <w3S(q,w)>                                                 [moment 3]
%     'm4'     M4   = <w4S(q,w)>                                                 [moment 4]
%     'sq_approx' S(q) ~ M2/Wl^2^2/sqrt(q)    structure factor estimate from pure inelastic
%
% output:
%   moments=[ sq M1 wc wl wq M2 M3 M4 ] as an iFunc model or array
%
% Example: m = moments(sqw_visco_elastic_simple,'wc'); plot(m);
% 

  if nargin < 2, select=[]; end
  if isempty(select), select = {'sq','recoil','wc','wl','wq','m2','m3','m4'}; end
  
  if iscellstr(select)
    M = [];
    for index=1:numel(select)
      M = [ M moments(self, select{index}) ];
    end
    return
  end
  
  % the resulting model has an energy axis
  M = copyobj(self);
  M.Dimension = 1;
  M.UserData.moment_type = select;
  select = lower(select);
  
  % check for classical model and temperature
  classical    = findfield(M, 'classical','first');
  is_classical = false;
  if ~isempty(classical)
    classical = get(self, classical);
    if classical(1), is_classical = true; end
  end
  if ~is_classical || any(strfind('wq normalized mean wc collective isothermal', select))
    gT = [ 'T = p(' num2str(numel(M.Parameters)+1) '); ' ];
    if isvector(M.Guess) && isnumeric(M.Guess)
      M.Guess = [ M.Guess(:)' 300 ];  % typical
    else
      if ~iscell(M.Guess), M.Guess = { M.Guess }; end
      M.Guess{end+1} = [ 300 ];
    end
    M.Parameters{end+1} = [ 'Temperature [K] ' mfilename ]; 
  else gT = [ 'T = 0; ' ];
  end

  % add Mass parameter
  if isvector(M.Guess) && isnumeric(M.Guess)
    M.Guess(end+1) = 12;
  else
    if ~iscell(M.Guess), M.Guess = { M.Guess }; end
    M.Guess{end+1} = [ 12];
  end
  M.Parameters{end+1} = [ 'Mass Material molar weight [g/mol] ' mfilename ];
  
  % we need a w axis. x axis is given as momentum when evaluating the moment.
  M = {[ 'q = x; q_sav_moment=q; ' gT ];
    ['kT=T/11.604; M=p(' num2str(numel(M.Parameters)) '); M1=0;'];
    'qmax=max(abs(q)); Emax=2.0721*qmax^2;';
    'w=linspace(-Emax, Emax, 200); this.UserData.moment_Emax = Emax;';
    '[q,w] = meshgrid(unique(q),w); x=q; y=w;'; } + M;
  % then evaluate the initial 2D model (from plus operator)
  
  % and append the moment computation
  M.Expression{end+1} = 'signal0=signal;';
  if any(strfind('sq M0 structure 0th moment wq normalized mean wc collective isothermal', select))
    M.Expression{end+1} = 'M0 = trapz(signal0); sq = M0;';
    M.Expression{end+1} = 'signal = sq;';
    M.Name=[ 'S(q) structure factor; ' self.Name ];
  end
  if any(strfind('m2 2nd moment wc collective isothermal wl harmonic longitudinal sq_approx', select))
    M.Expression{end+1} = 'M2 = trapz(w.^2.*signal0);';
    M.Expression{end+1} = 'signal = M2;';
    M.Name=[ '<w2S> 2nd moment; ' self.Name ];
  end
  if any(strfind('m4 4th moment', select))
    M.Expression{end+1} = 'M4 = trapz(w.^4.*signal0);';
    M.Expression{end+1} = 'signal = M4;';
    M.Name=[ '<w4S> 4th moment; ' self.Name ];
  end
  if any(strfind('m1 1st moment wr recoil wc collective isothermal wl harmonic longitudinal', select))
    M.Expression{end+1} = 'M1 = trapz(w.*signal0);';
    M.Expression{end+1} = 'signal = M1;';
    M.Name=[ '<wS> 1st moment recoil E_r=h^2q^2/2M; ' self.Name ];
  end
  if any(strfind('m3 3rd moment wl harmonic longitudinal', select))
    M.Expression{end+1} = 'M3 = trapz(w.^3.*signal0);';
    M.Expression{end+1} = 'signal = M3;';
    M.Name=[ '<w3S> 3rd moment; ' self.Name ];
  end
  if any(strfind('wq normalized mean', select))  % requires M0
    % half width from normalized 2nd frequency moment J-P.Hansen and I.R.McDonald 
    % Theory of simple liquids Academic Press New York 2006.
    % Lovesey p180 eq. 5.38 = w0
    M.Expression{end+1} = 'wq = 2*q.*sqrt(kT/M./M0);';
    M.Expression{end+1} = 'signal = wq;';
    M.Name=[ 'w_q=q sqrt[kB T/M S(q)] mean energy transfer; ' self.Name ];
  end
  
  if is_classical && any(strfind('wc collective isothermal', select))  % requires M0 M2
    % M2 = q.^2.*kT/M
    % sqrt(<w2S>/s(q)) == q sqrt(kT/M/s(q)) collective/isothermal
    M.Expression{end+1} = 'wc      = sqrt(M2./M0);';
    M.Expression{end+1} = 'signal  = wc;';
    M.Name=[ 'w_c collective/isothermal dispersion; ' self.Name ];
  end
  if is_classical && any(strfind('wl harmonic longitudinal', select))  % requires M2
    % maxima wL(q) of the longitudinal current correlation function ~ wl
    M.Expression{end+1} = 'm3      = abs(trapz(abs(w).^3.*signal0));';
    M.Expression{end+1} = 'wl      = m3./M2;';
    M.Expression{end+1} = 'signal  = wl;';
    M.Name=[ 'w_l harmonic/longitudinal excitation; ' self.Name ];
  end
  if ~is_classical && any(strfind('wc collective isothermal', select))  % requires M0 M1
    M.Expression{end+1} = 'wc      = sqrt(2*kT.*M1./M0);';
    M.Expression{end+1} = 'signal  = wc;';
    M.Name=[ 'w_c collective/isothermal dispersion; ' self.Name ];
  end
  if ~is_classical && any(strfind('wl harmonic longitudinal sq_approx', select))  % requires M3 M1
    M.Expression{end+1} = 'wl      = sqrt(M3./M1);';
    M.Expression{end+1} = 'signal  = wl;';
    M.Name=[ 'w_l harmonic/longitudinal excitation; ' self.Name ];
  end
  if any(strfind('sq_approx', select)) % requires wl M2
    % a very crude estimate of S(q) may be obtained from (phenomenological):
    %   M2./(wl-min(wl)).^2./sqrt(q) which works remarkably !
    % and normalize it.
    M.Expression{end+1} = 'sq = M2./(wl-min(wl)).^2./sqrt(q);';
    M.Expression{end+1} = 'signal = sq;';
    M.Name=[ 'S(q) structure factor estimate; ' self.Name ];
  end
  
  % restore initial axes
  M.Expression{end+1} = 'x = q_sav_moment';

  M = iFunc(M);
end
