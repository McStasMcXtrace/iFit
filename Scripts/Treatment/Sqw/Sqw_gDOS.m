function [gDOS, g] = Sqw_gDOS(s, method, n)
% Sqw_gDOS: compute the generalised density of states (gDOS)
%  The gDOS is an approximation of the vibrational spectra (DOS).
%  This routine should be applied on an incoherent dynamic S(q,w) data set.
%  The S(q,w) is a dynamic structure factor aka scattering function.
%
%       gDOS(q,w) = S(q,w) w/q2 [1 - exp(-hw/kT)] [Carpenter/Price]
%
%  The applicability to a coherent dynamic structure factor S(q,w) should be
%    taken with great care, as this formalism then does not hold.
%
%  The method to use in the gDOS computation can be given as 2nd argument
%       gDOS = Sqw_gDOS(Sqw, 'Carpenter')
%       gDOS = Sqw_gDOS(Sqw, 'Bellisent')
%       gDOS = Sqw_gDOS(Sqw, 'Bredov') better for coherent scatterers
%
%  The gDOS(w) is obtained by extracting the low momentum values out of gDOS(q,w).
%  The syntax is:
%       [g(w), g(q,w)]=Sqw_gDOS(Sqw, method, n)
%
% input:
%   s: Sqw data set e.g. 2D data set with w as 1st axis (rows, meV), q as 2nd axis (Angs-1).
%   method: 'Carpenter' (default),'Bellisent' or 'Bredov'
%   n: number of low-angle values to integrate (integer). Default is 10 when omitted.
%
% output:
%   g:   gDOS(w)   (1D iData versus energy)
%   g2D: gDOS(q,w) (2D iData)
%
% conventions:
% omega = Ei-Ef = energy lost by the neutron
%    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%    omega < 0, neutron gains energy, anti-Stokes
%
% references: Price J. et al, Non Cryst Sol 92 (1987) 153
%         Bellisent-Funel et al, J. Mol. Struct. 250 (1991) 213
%         Carpenter and Pelizarri, Phys. Rev. B 12, 2391 (1975)
%         Suck et al, Journal of Alloys and Compounds 342 (2002) 314
%         Bredov et al., Sov. Phys. Solid State 9, 214 (1967)
%
% Example: Sqw=iData('SQW_coh_lGe.nc'); g = Sqw_gDOS(Sqw_Bosify(Sqw_symmetrize(Sqw))); plot(g);

  g = []; gDOS=[];
  if nargin == 0, return; end
  if ischar(s)
    s = iData(s);
  end
  if ~isa(s, 'iData')
    disp([ mfilename ': ERROR: The data set should be an iData object, and not a ' class(s) ]);
    return; 
  end
  
  if nargin < 2, method = 'Carpenter'; end
  if nargin < 3, n=[]; end
  if isempty(n) || n <= 0, n=10; end
  
  % handle array of objects
  if numel(s) > 1
    for index=1:numel(s)
      [gDOS1, g1] = feval(mfilename, s(index), n);
      gDOS = [ gDOS gDOS1 ];
      g    = [ g g1 ];
    end
    return
  end
  
  s = Sqw_check(s);
  if isempty(s), return; end

  % test if classical
  if isfield(s,'classical') || ~isempty(findfield(s, 'classical'))
    if s.classical == 1
      disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source ' seems to be classical.' ])
      disp('  The gDOS computation may be wrong. Apply Sqw_Bosify first.');
    end
  end
  
  % compute g(q,w) aka P(alpha,beta) -------------------------------------------
  yl=getaxis(s, '1');
  xl=getaxis(s, '2');
  w= s{1}; 
  
  switch lower(method)
  case {'carpenter','price'}
    g = Sqw_gDOS_Carpenter(s);
  case 'bredov'
    g = Sqw_gDOS_Bredov(s);
  otherwise % 'bellisent'
    method = 'Bellisent';
    g = Sqw_gDOS_Bellisent(s);
  end

  if ischar(xl), setaxis(g,2, xl); end
  if ischar(yl), setaxis(g,1, yl); end

  title(g,[ 'gDOS 2D(' s.Title ') ' method ]);
  if isempty(g.Label), g.Label='gDOS'; end
  
  % this is the weighting for valid data
  g(~isfinite(g)) = 0;
  % normalise to 1. No need for the mass, temperature and other constant factors
  sum_g = g{0};
  g = g./sum(sum_g(:)); 
  
  % determine the low-momentum limit -------------------------------------------
  s0 = subsref(g,struct('type','.','subs','Signal'));
  gDOS = zeros(size(s0,1),1);  % column of g(w) DOS for q->0
  nc = min(n,size(s0,2)); % get the first n=10 values for each energy transfer
  for i=1:size(s0,1)
    nz = find(s0(i,:) > 0, nc, 'first');
    if ~isempty(nz)
      gDOS(i) = sum(s0(i,nz));
    end
  end
  gDOS=gDOS./sum(gDOS);
  if isvector(w)
    gDOS=iData(w,gDOS);
  else
    gDOS=iData(w(1,:),gDOS);
  end
  gDOS.Error=0;
  gDOS.Title = [ 'gDOS(' s.Title ') ' method ];
  title(gDOS,[ 'gDOS(' s.Title ') ' method ]);
  xlabel(gDOS, 'Energy');
  gDOS.Label=[ 'gDOS ' method ];






% ------------------------------------------------------------------------------
function g = Sqw_gDOS_Carpenter(s)
% See e.g.: Boatner et al, PRB 56 (1997) 11584
%           Price J. et al, Non Cryst Sol 92 (1987) 153
%           J. M. Carpenter and C. A. Pelizarri, Phys. Rev. B 12, 2391 (1975)
%
% g(q,w) = exp(2*W) ./ q.^2 * w./(n(w)+1) S(q,w)   [Carpenter]
%
% with:
%   exp(2*W) ~ exp(-0.003*q.^2) ~ 0.9-1 in q=0:5 Angs-1. Assumed to be 1.
%   n(hw) = 1./(exp(w/kT)-1))

  hw= s{1}; 
  q = s{2};
  T = Sqw_getT(s);
  if isempty(T)
    T = 293;
    disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source  ])
    disp(['    has no temperature defined. Using T=' num2str(T) ' K' ])
  end
  beta    = -11.605*hw/T;
  n = 1./(exp(beta) - 1);
  n(~isfinite(n)) = 0;
  
  g       = abs(hw./(n+1)).*s./q.^2;
  g.Error = s.Error./s.Signal.*g.Signal;
  
function g = Sqw_gDOS_Bellisent(s)
% See e.g.: Bellisent-Funel et al, J. Mol. Struct. 250 (1991) 213
%
% uses Carpenter with approx: exp(2W)=1 and 1/(n(w)+1) ~ hw
%
% g(q,w) ~ hw.^2./q.^2 S(q,w)                      [Bellisent]
  hw = s{1}; 
  q  = s{2};
  
  g = s.*hw.^2./q.^2;

function g = Sqw_gDOS_Bredov(s)
% See e.g.: Suck et al, Journal of Alloys and Compounds 342 (2002) 314
%           M. M. Bredov et al., Sov. Phys. Solid State 9, 214 (1967).
%
% g(w) = Ei w./(n(w)+1)/(q_max^4-q_min^4)
%        * \int_{theta_min -> theta_max} exp(2*W) kf/ki S(q,w) sin(theta) dtheta

  g = Sqw_gDOS_Carpenter(s);
  q = s{2};
  % for each energy, we get the min and max which can be measured
  % requires to know the incident energy.
  q4= zeros(size(s{1}));
  for index=1:size(s, 1)
    sw = s(index,:); % slab for a given energy. length is 'q'
    sw = sw{0};
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
