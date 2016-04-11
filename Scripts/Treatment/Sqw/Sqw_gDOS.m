function [gDOS, g] = Sqw_gDOS(s, n)
% Sqw_gDOS: compute the generalised density of states (gDOS)
%  The gDOS is an approximation of the vibrational spectra (DOS).
%  This routine should be applied on an incoherent dynamic S(q,w) data set.
%  The S(q,w) is a dynamic structure factor aka scattering function.
%
%               gDOS(q,w) = S(q,w) w/q2 [1 - exp(-hw/kT)]
%
%  The applicability to a coherent dynamic structure factor S(q,w) should be
%    taken with great care, as this formalism then does not hold.
%
%  The gDOS(w) is obtained by extracting the low momentum values out of gDOS(q,w).
%  The syntax is:
%                 [g(w), g(q,w)]=Sqw_gDOS(Sqw, n)
%
% input:
%   s: Sqw data set e.g. 2D data set with w as 1st axis (rows, meV), q as 2nd axis (Angs-1).
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
%             Bellisent-Funel et al, J. Mol. Struct. 250 (1991) 213
%
% Example: g = Sqw_gDOS(s); plot(g);

  g = []; gDOS=[];
  if nargin == 0, return; end
  if ~isa(s, 'iData')
    disp([ mfilename ': ERROR: The data set should be an iData object, and not a ' class(s) ]);
    return; 
  end
  
  if nargin == 1, n=[]; end
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
  
  w = s{1}; yl=getaxis(s, '1');
  q = s{2}; xl=getaxis(s, '2');

  g = s.*w.^2./q.^2;  % from Price 1987
  if ischar(xl), setaxis(g,2, xl); end
  if ischar(yl), setaxis(g,1, yl); end
  
  % reset axes
  % g{1}=w; g{2}=q;

  title(g,[ 'gDOS 2D(' s.Title ')' ]);
  if isempty(g.Label), g.Label='gDOS'; end
  
  % this is the weighting for valid data
  g(~isfinite(g)) = 0;
  
  % determine the low-momentum limit
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
  gDOS.Title = [ 'gDOS(' s.Title ')' ];
  title(gDOS,[ 'gDOS(' s.Title ')' ]);
  xlabel(gDOS, 'Energy');
  gDOS.Label='gDOS';
  setalias(gDOS,'Pab', Sqw_pab(s));

  
function Pab = Sqw_pab(s)
  % compute the P(q,hw)
  hw= s{1}; 
  q = s{2};
  y = s{0};
  yerr = get(s,'Error');
  T = Sqw_getT(s);
  if isempty(T)
    T = 293;
    disp([ mfilename ': WARNING: The data set ' s.Tag ' ' s.Title ' from ' s.Source  ])
    disp(['    has no temperature defined. Using T=' num2str(T) ' K' ])
  end
  alpha = 2.0721*q.^2/T;
  beta  = 11.605*hw/T;
  Pab     = beta.*(exp(beta) - 1).*s./alpha;

