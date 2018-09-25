function [g, fig] = dos(s, method, varargin)
% iData_Sqw2D: dos: compute the generalised density of states (gDOS) from a S(q,w)
%
%   g = dos(s, method, n, T, DW)
%
% compute: iData_Sqw2D -> generalised Density of States gDOS [p=1]
%
%  The returned generalised density of states corresponds with the 1-phonon term in the
%  the incoherent Gaussian approximation. This density of states is normalised to 1.
%
%       gDOS(q,w) = S(q,w) w^2/q^2                   Bellissent
%       gDOS(q,w) = S(q,w) w  /q^2/[1 + n(hw)]       Carpenter/Price
%  and:
%       gDOS(w)   = lim(q->0) [ gDOS(q,w) ]
%
%       gDOS(q,w) = w*q*S(q,w)*exp(2W(q))/[Qmax^4 - Qmin^4]/(1+n(w)) Bredov/Oskotskii
%       gDOS(w)   = trapz(g, 2)
%
%  The Bredov/Oskotskii methodology provides the best gDOS estimate, using the
%    whole data set.
%
%  LIMITATIONS/WARNINGS:
%  The incoherent approximation states that the gDOS from an incoherent S(q,w) is 
%    roughly equal to that obtained from a coherent S(q,w). However, the 
%    applicability to a coherent dynamic structure factor S(q,w) should be
%    taken with great care, as this formalism then does not fully hold.
%  This implementation is in principle exact for an isotropic monoatomic material,
%    e.g. a liquid, powder, or cubic crystal. 
%  This routine should better be applied on an incoherent dynamic S(q,w) data set.
%
%  The method to use in the gDOS computation can be given as 2nd argument
%       gDOS = dos(Sqw, 'Bredov')         more accurate as it uses 100% of data
%       gDOS = dos(Sqw, 'Carpenter')      Temperature must be a property
%       gDOS = dos(Sqw, 'Bellissent')     simple yet efficient
%
%  The gDOS is stored in the 'gDOS' property of the initial Sqw object, and is 
%    retrieved without recomputation when available. 
%  To force a recomputation of the gDOS, use:
%       dos(Sqw, 'method', 0) or dos(Sqw, 'method', 'force')
%
%  Input arguments can be given in order, or with name-value pairs, or as a 
%    structure with named fields.
%
% References: Price J. et al, Non Cryst Sol 92 (1987) 153
%         Bellissent-Funel et al, J. Mol. Struct. 250 (1991) 213
%         Carpenter and Pelizarri, Phys. Rev. B 12, 2391 (1975)
%         Suck et al, Journal of Alloys and Compounds 342 (2002) 314
%         Bredov et al., Sov. Phys. Solid State 9, 214 (1967)
%         V.S. Oskotskii, Sov. Phys. Solid State 9 (1967), 420.
%         H. Schober, Journal of Neutron Research 17 (2014) 109–357. (esp. p307-315)
%
% syntax:
%   g = dos(s)
%   g = dos(s, method, n, T, DW)
%   g = dos(s, 'method', method, 'n', n, 'T', T, 'DW', dw)
%
% input:
%   s:      Sqw data set e.g. 2D data set with w as 1st axis (rows, meV), q as 2nd axis (Angs-1).
%   method: 'Carpenter','Bellissent' or 'Bredov' (default)
%   n:      number of low-angle values to integrate (integer). Default is 10 when omitted.
%           when 0 or 'force', the gDOS is re-computed. Only for Carpenter/Bellissent.
%   T:      optional temperature to use for computation (leave undefined for automatic).
%   DW:     optional Debye-Waller coefficient gamma=<u^2> [Angs^2] e.g. 0.005
%           The Debye-Waller function is      2W(q)=gamma*q^2
%           The Debye-Waller factor   is exp(-2W(q))
%
% output:
%   g:      gDOS(w)   (1D iData versus energy)
%
% Example: Sqw=iData_Sqw2D('D2O_liq_290_coh.sqw.zip'); g = dos(Bosify(symmetrize(Sqw))); plot(g);
%
% See also: iData_Sqw2D/multi_phonons, iData_Sqw2D/incoherent
%           iData_vDOS/multi_phonons, iData_vDOS/multi_phonons
% (c) E.Farhi, ILL. License: EUPL.

  g=[]; fig = [];
  if isempty(s), return; end
  
  if nargin < 2, method = '';
  elseif ~any(strcmpi(method, {'force','bredov','oskotskii','default','carpenter','bellissent'})) && ~isempty(method)
    varargin = { method varargin{:} };
    method = '';
  end
  if isempty(method), method = 'default'; end
  
  p = varargin2struct({'n' 't' 'DW' 'gamma' 'u2'}, varargin, true);
  if ~isfield(p,'method') p.method = method; end
  
  if isfield(p, 'gamma') && ~isempty(p.gamma), p.dw = p.gamma; end
  if isfield(p, 'u2')    && ~isempty(p.u2),    p.dw = p.u2; end

  % handle array of objects
  if numel(s) > 1
    for index=1:numel(s)
      g = [ g feval(mfilename, s(index), p) ];
    end
    if ~isempty(inputname(1))
      assignin('caller',inputname(1),s);
    end
    return
  end

  % check for method arguments
  if isempty(s), return; end
  if strcmp(p.method, 'force'), p.method = []; p.n=0; end
  if isempty(p.method),         p.method = 'default'; end
  if ~isempty(p.n) && (strcmp(p.n, 'force') || (isnumeric(p.n) && p.n <= 0))
    s = rmalias(s, 'gDOS');
    p.n = [];
  end
  
  if isfield(s, 'gDOS') && nargin == 1
    % access the stored result (cached)
    g = get(s, 'gDOS');
  else
    % do the computation
    
    if isempty(p.t)
      p.t = Sqw_getT(s);
    end
    if isempty(p.t) && ~isempty(strfind(p.method, 'Carpenter')) 
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source  ])
      disp(['    has no Temperature defined. Will assume 1/[1+p.n(w)] ~ w.' ])
    end
    
    switch lower(p.method)
    case {'bredov','oskotskii','default'}
      [g,w] = sqw_phonon_dos_Bredov(s, p.t, p.dw); 
      method = 'Bredov';
    otherwise % {'carpenter','price','bellissent'}
      if isempty(p.n), p.n=max(10, size(s, 2)/10); end
      [g, p.method] = sqw_phonon_dos_Carpenter(s, p.t, p.n, p.dw); % also incl. Bellissent
    end

    % check if DOS is on both energy axes. If so, symmetrize.
    w = getaxis(g,1);
    if any(w(:) < 0) && any(w(:) > 0)
      % fold DOS on positive energy side
      g = combine(g, setaxis(g, 1, -getaxis(g,1)));
      g = xlim(g, [0 inf]);
    elseif any(w(:) < 0)
      % set energy axis positive
      g = setaxis(g, 1, -getaxis(g,1));
    end
    
    
    % get valid data
    index = find(~isfinite(g));
    g0    = getaxis(g, 0); g0(index) = 0;
    g     = subsasgn(g, struct('type','.','subs','Signal'), g0);

    % make it a vDOS object and normalise
    g     = iData_vDOS(g);
    g     = g/sum(g);
    
    % set Title, etc
    g.Title = [ 'DOS(' s.Title ') ' p.method ];
    g.Label = [ 'gDOS ' p.method ];
    title( g,[ 'gDOS [p=1](' s.Title ') ' p.method ]);
    xlabel(g,  ylabel(s));  % energy axis
    g = commandhistory(g,    'dos', s, p.method, p.n);
    setalias(g, 'DOS_method', p.method, 'Method used for DOS estimate');
    setalias(s, 'gDOS', g, g.Title);
    
    if length(inputname(1))
      assignin('caller',inputname(1),s);
    end
  end % compute

  % plot total DOS when no output
  if nargout == 0 && ~isempty(g)
    fig=figure; 
    h=plot(g); set(h,'LineWidth',2);
    set(fig, 'NextPlot','new');
  end

% ------------------------------------------------------------------------------
function [g, method] = sqw_phonon_dos_Carpenter(s, T, nc, DW)
% See e.g.: Boatner et al, PRB 56 (1997) 11584
%           Price J. et al, Non Cryst Sol 92 (1987) 153
%           J. M. Carpenter and C. A. Pelizarri, Phys. Rev. B 12, 2391 (1975)
% See e.g.: Bellissent-Funel et al, J. Mol. Struct. 250 (1991) 213
%
% g(q,w) = exp(2*W) ./ q.^2 * w./(n(w)+1) S(q,w)   [Carpenter]
% g(q,w) ~ hw.^2./q.^2 S(q,w)                      [Bellissent]
%
% with:
%   exp(2*W) ~ exp(-0.003*q.^2) ~ 0.9-1 in q=0:5 Angs-1. Assumed to be 1.
%   n(hw) = 1./(exp(w/kT)-1))

  if isfield(s,'classical') || ~isempty(findfield(s, 'classical'))
    classical = get(s,'classical');
  else classical = [];
  end
  
  hw= getaxis(s,1); 
  q = getaxis(s,2);
  if ~isempty(DW)
    DW= exp(2*DW*q.^2);
  else DW = 1;
  end
  
  % compute g(q,w) aka P(alpha,beta) 
  if (~isempty(classical) && classical(1)) || isempty(T) || T==0 % do not use Temperature
    % Bellissent
    % uses Carpenter with approx: exp(2W)=1 and 1/(p.n(w)+1) ~ hw
    method = 'Bellissent';
    
    g = s.*hw.^2./(q.^2);
  else
    % Carpenter/Price
    method = 'Carpenter';

    beta    = 11.605*hw/T;
    n = 1./(exp(beta) - 1);
    n(~isfinite(n)) = 0;
    g       = abs(hw./(n+1)).*s./(q.^2);
  end
  g = g.*DW;

  % determine the low-momentum limit -------------------------------------------
  s0  = subsref(g,struct('type','.','subs','Signal'));
  DOS = zeros(size(s0,1),1);  % column of g(w) DOS for q->0
  
  nc  = ceil(min(nc,size(s0,2)));    % get the first nc values for each energy transfer
  for i=1:size(s0,1)
    nz = find(isfinite(s0(i,:)) & s0(i,:) > 0, nc, 'first');
    if ~isempty(nz)
      DOS(i) = sum(s0(i,nz));
    end
  end

  % set new gDOS value
  rmaxis( g, 2);
  g = subsasgn(g, struct('type','.','subs','Signal'), DOS);
  setaxis(g, 1, unique(hw));

% ------------------------------------------------------------------------------
function [g, w] = sqw_phonon_dos_Bredov(s, T, DW)
% See e.g.: Suck et al, Journal of Alloys and Compounds 342 (2002) 314
%           M. M. Bredov et al., Sov. Phys. Solid State 9, 214 (1967).
%
% g(w) = Ei w./(n(w)+1)/(q_max^4-q_min^4) m/sigma exp(2*W) 
%        * \int_{q_min -> q_max} q / kf DDCS dq
%
%      = w./(n(w)+1)/(q_max^4-q_min^4) m/sigma exp(2*W) * \int_{q_min -> q_max} q S(q,w) dq

% Schober (9.299) p 315. See also (10.96)
%   \int sin(theta) DDCS d(theta) = sigma ki^2/m exp(-2W(q)) (q^4max-q^4min) g(w) (n+1)/w
%   \int q/ki/kf    DDCS dq
% and (usual Squires exp)
%   DDCS(q,w) = d2sigma/dOmega/dEf = kf/ki sigma S(q,w)

%   \int q/ki/kf kf/ki sigma S(q,w) dq = sigma ki^2/m exp(-2W(q)) (q^4max-q^4min) g(w) (n+1)/w
%   \int q/ki^2 S(q,w) dq              =       ki^2/m exp(-2W(q)) (q^4max-q^4min) g(w) (n+1)/w
% with fixed Ki:
%   \int q S(q,w) dq                   = [ exp(-2W(q))/m (q^4max-q^4min)  (p.n+1)/w ] g(w)
%
% g(w) =  \int q S(q,w) dq / [ exp(-2W(q))/m (q^4max-q^4min)  (n+1)/w ]
%      = [\int q S(q,w) dq] exp(2W(q)) m / (q^4max-q^4min) * w/(1+n)
%      ~ [\int q S(q,w) dq] exp(2W(q)) m w^2/(q^4max-q^4min)

% w/(1+n) ~ w2 to avoid divergence

  if isfield(s,'classical') || ~isempty(findfield(s, 'classical'))
    classical = get(s,'classical');
  else classical = [];
  end

  % get an estimate of corresponding incident energy
  [s,lambda,distance,chwidth,Ei,Ki] = Sqw_search_lambda(s);
  parameters = get(s, 'parameters');
  
  % restrict s to a dynamic range (so that q4 corresponds with a 'simulated' experiment)
  s = dynamic_range(s, Ei);
  
  % re-sample histogram when axes are not vectors
  hist_me = false;
  for index=1:ndims(s)
    x = getaxis(s, index);
    if numel(x) ~= length(x), hist_me=true; break; end
  end
  if hist_me
    disp([ mfilename ': INFO: Re-sampling on regular grid ' s.Tag ' ' s.Title ' from ' s.Source ])
    s   = hist(s, size(s));  % the data set is well covered here and interpolation works well.
  end
  
  % axes
  w = getaxis(s,1);
  q = getaxis(s,2);

  if ~isempty(DW)
    DW= exp(2*DW*q.^2);
  else DW = 1;
  end
  
  qSq = s.*q.*DW;
  
  
  % compute delta(Q)
  q4 = zeros(size(w));
  sd = double(s);
  qmax = []; qmin = [];
  for index=1:size(s, 1)
    sw = sd(index,:); % slab for a given energy. length is 'q'
    valid = find(isfinite(sw) & sw > 0);
    if isvector(q)
      qw = q;
    else
      qw = q(index,:);
    end
    q_min = min(qw(valid)); q_max = max(qw(valid));
    if isempty(q_max) || isempty(q_min)
      q_max = inf; q_min = 0;
    end
    qmax(end+1) = q_max;
    qmin(end+1) = q_min;
  end
  qmax = qmax(:); qmin = qmin(:);
  index = find(abs(qmax - qmin) < .1);
  qmax(index) = qmin(index) + 0.1;

  % compute the dos
  g = [];
  if (~isempty(classical) && classical(1)) || isempty(T) || T==0 % do not use Temperature
    % g(w) = [\int q S(q,w) dq] exp(2W(q)) m / (q^4max-q^4min) * w.^2
    g = trapz(qSq,2).*abs(w.^2./(qmax.^4 - qmin.^4));
  else    % use Bose/Temperature (in quantum case)
    beta    = 11.605*w/T;
    p.n = 1./(exp(beta) - 1);
    p.n(~isfinite(p.n)) = 0;
    % g(w) = [\int q S(q,w) dq] exp(2W(q)) m / (q^4max-q^4min) * w/(1+n)
    g = trapz(qSq,2).*abs(w./(p.n+1)./(qmax.^4 - qmin.^4));
  end
  
  % set back aliases
  setalias(g, 'parameters', parameters, 'Material parameters'); 

