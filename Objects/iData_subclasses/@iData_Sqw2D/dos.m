function [DOS, fig] = dos(s, method, n)
% iData_Sqw2D: dos: compute the generalised density of states (gDOS)
%
%  The gDOS is an approximation of the vibrational spectra (DOS).
%  This routine should better be applied on an incoherent dynamic S(q,w) data set.
%  The incoherent approximation states that the gDOS from an incoherent S(q,w) is 
%  roughly equal to that obtained from a coherent S(q,w).
%
%       gDOS(q,w) = S(q,w) w^2/q^2                   Bellissent
%       gDOS(q,w) = S(q,w) w  /q^2/[1 + n(hw)]       Carpenter/Price
%  and:
%       gDOS(w)   = lim(q->0) [ gDOS(q,w) ]
%
%  The applicability to a coherent dynamic structure factor S(q,w) should be
%  taken with great care, as this formalism then does not fully hold.
%
%  The method to use in the gDOS computation can be given as 2nd argument
%       gDOS = dos(Sqw, 'Carpenter')      Temperature must be a property
%       gDOS = dos(Sqw, 'Bellissent')     simple yet efficient
%       gDOS = dos(Sqw, 'Bredov')         better for coherent S(q,w)
%
%  The gDOS is stored in the 'gDOS' property, and is retrieved without recomputation
%  when available. To force a recomputation of the gDOS, use:
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
%   method: 'Carpenter','Bellissent' (default) or 'Bredov'
%   n:      number of low-angle values to integrate (integer). Default is 10 when omitted.
%           when 0 or 'force', the gDOS is re-computed.
%
% output:
%   g:      gDOS(w)   (1D iData versus energy)
%
% Example: Sqw=iData_Sqw2D('SQW_coh_lGe.nc'); g = dos(Bosify(symmetrize(Sqw))); plot(g);
%
% See also: iData_Sqw2D/multi_phonons, iData_Sqw2D/incoherent
%           iData_vDOS/multi_phonons_dos, iData_vDOS/multi_phonons_incoherent
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
  if isempty(method), method = 'Bellissent'; end
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
      if get(s,'classical')
        disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' seems to be classical.' ])
        disp('  The gDOS computation may be wrong.');
      end
    end
    
    % compute g(q,w) aka P(alpha,beta) -------------------------------------------
    yl=getaxis(s, '1');
    xl=getaxis(s, '2');
    w= getaxis(s,1);
    
    T = Sqw_getT(s);
    if isempty(T) && isempty(strfind(method, 'Bellissent'))
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source  ])
      disp(['    has no Temperature defined. Using method=''Bellissent''' ])
      method = 'Bellissent';
    end
    
    switch lower(method)
    case {'carpenter','price'}          
      g = sqw_phonon_dos_Carpenter(s, T);  % requires Temperature 
    case 'bredov'
      g = sqw_phonon_dos_Bredov(s, T);     % requires Temperature (call Carpenter)
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
    nc  = min(n,size(s0,2)); % get the first n=10 values for each energy transfer
    for i=1:size(s0,1)
      nz = find(s0(i,:) > 0, nc, 'first');
      if ~isempty(nz)
        DOS(i) = sum(s0(i,nz));
      end
    end
    DOS=DOS./sum(DOS);
    if isvector(w)
      DOS=iData(w,DOS);
    else
      DOS=iData(w(1,:),DOS);
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
    
    DOS       = iData_vDOS(DOS);
  end
  
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

function g = sqw_phonon_dos_Bredov(s, T)
% See e.g.: Suck et al, Journal of Alloys and Compounds 342 (2002) 314
%           M. M. Bredov et al., Sov. Phys. Solid State 9, 214 (1967).
%
% g(w) = Ei w./(n(w)+1)/(q_max^4-q_min^4)
%        * \int_{theta_min -> theta_max} exp(2*W) kf/ki S(q,w) sin(theta) dtheta

  g = sqw_phonon_dos_Carpenter(s,T);
  q = getaxis(s,2); 
  % for each energy, we get the min and max which can be measured
  % requires to know the incident energy.
  q4= zeros(size(getaxis(s,1)));
  s = double(s);
  for index=1:size(s, 1)
    sw = s(index,:); % slab for a given energy. length is 'q'
    valid = find(sw(isfinite(sw) & sw > 0));
    if isvector(q)
      qw = q;
    else
      qw = q(index,:);
    end
    q_min = min(qw(valid)); q_max = max(qw(valid));
    if isvector(q)
      q4(index) = q_max^4 - q_min^4;
    else
      q4(index,:) = q_max^4 - q_min^4;
    end
  end
  g       = g.*q.^2./q4;
  
