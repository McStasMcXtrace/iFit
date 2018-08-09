function [DOS, g, fig] = dos(s, method, n)
% iData_Sqw2D: dos: compute the generalised density of states (gDOS) from a S(q,w)
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
%       gDOS(q,w) = w*q*S(q,w)/[Qmax^4 - Qmin^4]/(1+n(w)) Bredov/Oskotskii
%       gDOS(w)   = trapz(g, 2)
%
%  The Bredov/Oskotskii methodology provides the best gDOS estimate, using the
%  whole data set.
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
% References: Price J. et al, Non Cryst Sol 92 (1987) 153
%         Bellissent-Funel et al, J. Mol. Struct. 250 (1991) 213
%         Carpenter and Pelizarri, Phys. Rev. B 12, 2391 (1975)
%         Suck et al, Journal of Alloys and Compounds 342 (2002) 314
%         Bredov et al., Sov. Phys. Solid State 9, 214 (1967)
%         V.S. Oskotskii, Sov. Phys. Solid State 9 (1967), 420.
%         H. Schober, Journal of Neutron Research 17 (2014) 109–357. (esp. p307)
%
% syntax:
%   g = dos(s)
%   g = dos(s, method, n)
%
% input:
%   s:      Sqw data set e.g. 2D data set with w as 1st axis (rows, meV), q as 2nd axis (Angs-1).
%   method: 'Carpenter','Bellissent' or 'Bredov' (default)
%   n:      number of low-angle values to integrate (integer). Default is 10 when omitted.
%           when 0 or 'force', the gDOS is re-computed.
%
% output:
%   g:      gDOS(w)   (1D iData versus energy)
%
% Example: Sqw=iData_Sqw2D('SQW_coh_lGe.nc'); g = dos(Bosify(symmetrize(Sqw))); plot(g);
%
% See also: iData_Sqw2D/multi_phonons, iData_Sqw2D/incoherent
%           iData_vDOS/multi_phonons, iData_vDOS/multi_phonons
% (c) E.Farhi, ILL. License: EUPL.

  DOS=[]; fig = [];
  if nargin < 2, method = []; end
  if nargin < 3, n=[]; end
  
  % handle array of objects
  if numel(s) > 1
    for index=1:numel(s)
      DOS = [ DOS feval(mfilename, s(index), method, n) ];
    end
    if ~isempty(inputname(1))
      assignin('caller',inputname(1),s);
    end
    return
  end

  if strcmp(method, 'force'), method = []; n=0; end
  if isempty(method), method = 'Bredov'; end
  if ~isempty(n) && (strcmp(n, 'force') || (isnumeric(n) && n <= 0))
    s = rmalias(s, 'gDOS');
    n = [];
  end
  if isempty(n), n=10; end
  if isempty(s), return; end
  
  if isfield(s, 'gDOS')
    % access the stored result 9cached)
    DOS = get(s, 'gDOS');
  else
    % do the computation
    % test if classical
    if isfield(s,'classical') || ~isempty(findfield(s, 'classical'))
      classical = get(s,'classical');
      if classical(1) && strcmp(lower(method), 'carpenter')
        method = 'Bellissent';  % do not use Bose/temp on classical
      end
    end
    
    % compute g(q,w) aka P(alpha,beta) -------------------------------------------
    yl=getaxis(s, '1');
    xl=getaxis(s, '2');
    w= getaxis(s,1);
    
    T = Sqw_getT(s);
    if isempty(T) && isempty(strfind(method, 'Bellissent')) && isempty(strfind(method, 'Bredov'))
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source  ])
      disp(['    has no Temperature defined. Using method=''Bellissent'' as it does not require it.' ])
      method = 'Bellissent';
    end
    
    switch lower(method)
    case {'carpenter','price'}          
      g = sqw_phonon_dos_Carpenter(s, T);  % requires Temperature 
    case {'bredov','oskotskii'}
      [g,w] = sqw_phonon_dos_Bredov(s, T);     % requires Temperature
      n = inf;                             % we use all the data.
    otherwise % 'bellissent'
      method = 'Bellissent';
      g = sqw_phonon_dos_Bellissent(s);
    end

    if ischar(xl), setaxis(g,2, xl); end
    if ischar(yl), setaxis(g,1, yl); end

    title(g,[ 'gDOS 2D(' s.Title ') ' method ]);
    if isempty(g.Label), g.Label='gDOS'; end
    
    % this is the weighting for valid data
    S.type='()';
    S.subs={ ~isfinite(g) };
    g = subsasgn(g, S, 0);
    % normalise to 1. No need for the mass, temperature and other constant factors
    sum_g = getaxis(g,0);
    g = g./sum(sum_g(:)); 
    
    % determine the low-momentum limit -------------------------------------------
    s0  = subsref(g,struct('type','.','subs','Signal'));
    DOS = zeros(size(s0,1),1);  % column of g(w) DOS for q->0
    
    if ~isfinite(n) || n<=0       % Bredov case
      DOS = sum(s0,2);
    else
      nc  = min(n,size(s0,2));    % get the first n=10 values for each energy transfer
      for i=1:size(s0,1)
        nz = find(isfinite(s0(i,:)) & s0(i,:) > 0, nc, 'first');
        if ~isempty(nz)
          DOS(i) = sum(s0(i,nz));
        end
      end
    end
    DOS=DOS./sum(DOS);        % normalise to 1
    if isvector(w)
      DOS=iData(w,DOS);
    else
      DOS=iData(w(1,:),DOS);  % from meshgrid
    end
    set(DOS,'Error',0);
    DOS.Title = [ 'DOS(' s.Title ') ' method ];
    DOS.Label = [ 'gDOS ' method ];
    title(DOS,[ 'gDOS(' s.Title ') ' method ]);
    xlabel(DOS, 'Energy');
    
    % copy initial aliases and UserData
    DOS.UserData = s.UserData;
    f = getalias(s);
    for index=1:numel(getalias(s))
      if ~isfield(DOS, f{index})
        [link, lab] = getalias(s, f{index});
        DOS = setalias(DOS, f{index}, link, lab);
      end
    end
    
    if numel(getaxis(DOS,1)) ~= prod(size(DOS))
      w = unique(getaxis(s, 1));
      DOS = setaxis(DOS, 1, w);
    end
    
    DOS       = iData_vDOS(DOS);
  end
  
  label(DOS, 0, [  'gDOS [p=1]' '(' label(s, 0) ')' ]);
  DOS = commandhistory(DOS,    'dos', s, method, n);
  
  if nargout == 0 && ~isempty(DOS)
    fig=figure; 
    % plot total DOS
    h=plot(DOS); set(h,'LineWidth',2);
    set(fig, 'NextPlot','new');
  end
  
  setalias(s, 'gDOS', DOS, DOS.Title);
  
  if nargout == 0 & length(inputname(1))
    assignin('caller',inputname(1),s);
  end
  

% ------------------------------------------------------------------------------
function g = sqw_phonon_dos_Carpenter(s, T)
% See e.g.: Boatner et al, PRB 56 (1997) 11584
%           Price J. et al, Non Cryst Sol 92 (1987) 153
%           J. M. Carpenter and C. A. Pelizarri, Phys. Rev. B 12, 2391 (1975)
%
% g(q,w) = exp(2*W) ./ q.^2 * w./(n(w)+1) S(q,w)   [Carpenter]
%
% with:
%   exp(2*W) ~ exp(-0.003*q.^2) ~ 0.9-1 in q=0:5 Angs-1. Assumed to be 1.
%   n(hw) = 1./(exp(w/kT)-1))

  hw= getaxis(s,1); 
  q = getaxis(s,2); 
  beta    = 11.605*hw/T;
  n = 1./(exp(beta) - 1);
  n(~isfinite(n)) = 0;
  
  g       = abs(hw./(n+1)).*s./q.^2;
  set(g,'Error', get(s,'Error')./get(s,'Signal').*get(g,'Signal'));
  
% ------------------------------------------------------------------------------
function g = sqw_phonon_dos_Bellissent(s)
% See e.g.: Bellissent-Funel et al, J. Mol. Struct. 250 (1991) 213
%
% uses Carpenter with approx: exp(2W)=1 and 1/(n(w)+1) ~ hw
%
% g(q,w) ~ hw.^2./q.^2 S(q,w)                      [Bellissent]
  hw= getaxis(s,1); 
  q = getaxis(s,2); 
  
  g = s.*hw.^2./q.^2;
  set(g,'Error', get(s,'Error')./get(s,'Signal').*get(g,'Signal'));

% ------------------------------------------------------------------------------
function [g, w] = sqw_phonon_dos_Bredov(s, T)
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
%   \int q S(q,w) dq                   = [ exp(-2W(q))/m (q^4max-q^4min)  (n+1)/w ] g(w)
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
  
  % restrict s to a dynamic range (so that q4 corresponds with a 'simulated' experiment)
  s = dynamic_range(s, Ei);
  
  % re-sample histogram when axes are not vectors
  hist_me = false;
  for index=1:ndims(s)
    x = getaxis(s, index);
    if numel(x) ~= length(x), hist_me=true; break; end
  end
  if hist_me
    s   = hist(s, size(s));  % the data set is well covered here and interpolation works well.
  end
  
  % axes
  w = getaxis(s,1);
  q = getaxis(s,2);
  
  qSq = q.*s;
  
  
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
  if ~isempty(classical) && classical(1) % do not use Temperature
    % g(w) = [\int q S(q,w) dq] exp(2W(q)) m / (q^4max-q^4min) * w.^2
    g = trapz(qSq,2).*abs(w.^2./(qmax.^4 - qmin.^4));
  else    % use Bose/Temperature (in quantum case)
    beta    = 11.605*w/T;
    n = 1./(exp(beta) - 1);
    n(~isfinite(n)) = 0;
    % g(w) = [\int q S(q,w) dq] exp(2W(q)) m / (q^4max-q^4min) * w/(1+n)
    g = trapz(qSq,2).*abs(w./(n+1)./(qmax.^4 - qmin.^4));
  end

