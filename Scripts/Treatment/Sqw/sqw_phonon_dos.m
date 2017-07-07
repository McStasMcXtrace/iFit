function [DOS, DOS_partials] = sqw_phonon_dos(s, method, n)
% sqw_phonon_dos: compute the density of states (gDOS or vDOS)
%
%  The routine can be used for 2D and 4D models and data sets.
%    when used on 4D data sets and models S(HKL,w), the vDOS is computed.
%    when used on 2D data sets and models S(|q|,w), the gDOS is computed.
%
%  ================================ 4D case ====================================
%    DOS = sqw_phonon_dos(s)    returns the vibrational density of states (vDOS)
%      the vDOS and the partials per mode are also stored in the UserData.
%    DOS = sqw_phonon_dos(s, n) does the same with n-bins on the vDOS (n=100)
%    when the DOS has already been computed, it is used as is. To force a
%    recomputation, set:
%      s.UserData.DOS=[];
%    to smooth the resulting distribution, use:
%      sDOS = smooth(DOS); plot(sDOS);
%  
%  ================================ 2D case ====================================
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
%       gDOS = sqw_phonon_dos(Sqw, 'Carpenter')
%       gDOS = sqw_phonon_dos(Sqw, 'Bellisent')
%       gDOS = sqw_phonon_dos(Sqw, 'Bredov') better for coherent scatterers
%
%  The gDOS(w) is obtained by extracting the low momentum values out of gDOS(q,w).
%  The syntax is:
%       [g(w), g(q,w)]=sqw_phonon_dos(Sqw, method, n)
%
% input:
%   s: Sqw data set e.g. 2D data set with w as 1st axis (rows, meV), q as 2nd axis (Angs-1).
%   method: 'Carpenter' (default),'Bellisent' or 'Bredov'
%   n: number of low-angle values to integrate (integer). Default is 10 when omitted.
%
% output:
%   g:   gDOS(w)   (1D iData versus energy)
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
% Example: Sqw=iData('SQW_coh_lGe.nc'); g = sqw_phonon_dos(Sqw_Bosify(Sqw_symmetrize(Sqw))); plot(g);
% (c) E.Farhi, ILL. License: EUPL.

  DOS=[]; DOS_partials=[];
  if nargin == 0, return; end
  
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
  
  if (isa(s,'iFunc') || isa(s,'iData')) && ndims(s) == 4
    [DOS, DOS_partials, gDOS, s] = sqw_phonon_dos_4D(s, method);
    if ~isempty(inputname(1))
      assignin('caller',inputname(1),s);
    end
    if nargout == 0 % plot
      fig=figure;
      DOS = s.UserData.DOS;
      DOS{1} = DOS{1}; % change energy unit
      xlabel(DOS,[ 'Energy [meV]' ]);
      % plot any partials first
      if isfield(s.UserData,'DOS_partials') && numel(s.UserData.DOS_partials) > 0
        d=s.UserData.DOS_partials;
        for index=1:numel(d)
          this_pDOS=d(index);
          this_pDOS{1} = this_pDOS{1};
          d(index) = this_pDOS;
        end
        h=plot(d);
        if iscell(h), h=cell2mat(h); end
        set(h,'LineStyle','--');
        hold on
      end
      % plot total DOS and rotate
      h=plot(DOS); set(h,'LineWidth',2);
    end
    return
  end
  
  if isempty(method), method = 'Carpenter'; end
  if isempty(n) || n <= 0, n=10; end

  if ~isa(s, 'iData'), s=iData(s); end

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
    g = sqw_phonon_dos_Carpenter(s);
  case 'bredov'
    g = sqw_phonon_dos_Bredov(s);
  otherwise % 'bellisent'
    method = 'Bellisent';
    g = sqw_phonon_dos_Bellisent(s);
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
  DOS = zeros(size(s0,1),1);  % column of g(w) DOS for q->0
  nc = min(n,size(s0,2)); % get the first n=10 values for each energy transfer
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
  DOS.Error=0;
  DOS.Title = [ 'DOS(' s.Title ') ' method ];
  title(DOS,[ 'gDOS(' s.Title ') ' method ]);
  xlabel(DOS, 'Energy');
  DOS.Label=[ 'gDOS ' method ];

% ------------------------------------------------------------------------------
function g = sqw_phonon_dos_Carpenter(s)
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
  
function g = sqw_phonon_dos_Bellisent(s)
% See e.g.: Bellisent-Funel et al, J. Mol. Struct. 250 (1991) 213
%
% uses Carpenter with approx: exp(2W)=1 and 1/(n(w)+1) ~ hw
%
% g(q,w) ~ hw.^2./q.^2 S(q,w)                      [Bellisent]
  hw = s{1}; 
  q  = s{2};
  
  g = s.*hw.^2./q.^2;

function g = sqw_phonon_dos_Bredov(s)
% See e.g.: Suck et al, Journal of Alloys and Compounds 342 (2002) 314
%           M. M. Bredov et al., Sov. Phys. Solid State 9, 214 (1967).
%
% g(w) = Ei w./(n(w)+1)/(q_max^4-q_min^4)
%        * \int_{theta_min -> theta_max} exp(2*W) kf/ki S(q,w) sin(theta) dtheta

  g = sqw_phonon_dos_Carpenter(s);
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
  
% ------------------------------------------------------------------------------

function [DOS, DOS_partials, gDOS, s] = sqw_phonon_dos_4D(s, n)
  % sqw_phonon_dos_4D: compute the phonon/vibrational density of states (vDOS)
  %
  % input:
  %   s: phonon S(q,w) (4D) Model or Data set (iFunc/iData)
  %
  % output:
  %   DOS:          phonon/vibrational density of states (iData)
  %   DOS_partials: phonon/vibrational density of state partials (iData array)
  %   gDOS:         integrated density of state, i.e. neutron weighted when 
  %                   extracted from a Model evaluation (iData)
  %
  % Example: DOS=sqw_phonon_dos(sqw_cubic_monoatomic('defaults'))
  % (c) E.Farhi, ILL. License: EUPL.

  DOS = []; gDOS=[]; f=[]; DOS_partials = [];
  
  % must be 4D iFunc or iData
  if ~nargin || (~isa(s, 'iFunc') && ~isa(s, 'iData')) || ndims(s) ~= 4
    return
  end
  
  if isa(s, 'iData')
    f = s;
  end
  if nargin < 2, n=[]; end
  
  % first get a quick estimate of the max frequency
  if  (~isfield(s.UserData,'DOS')  || isempty(s.UserData.DOS) ...
    || ~isfield(s.UserData,'gDOS') || isempty(s.UserData.gDOS)) && isa(s, 'iFunc')
    qh=linspace(-.5,.5,10);qk=qh; ql=qh; w=linspace(0.01,50,11);
    f=iData(s,[],qh,qk,ql',w);
    if isfield(s.UserData, 'FREQ') && ~isempty(s.UserData.FREQ)
      s.UserData.maxFreq = max(s.UserData.FREQ(:));
      disp([ mfilename ': maximum phonon energy ' num2str(max(s.UserData.maxFreq)) ' [meV] in ' s.Name ]);
    end
    if ~isfield(s.UserData, 'maxFreq') || isempty(s.UserData.maxFreq) ...
      || ~isfinite(s.UserData.maxFreq) || s.UserData.maxFreq <= 0
      s.UserData.maxFreq = 100;
    end
    
    % evaluate the 4D model onto a mesh filling the Brillouin zone [-0.5:0.5 ]
    s.UserData.DOS     = [];  % make sure we re-evaluate again on a finer grid
    s.UserData.maxFreq = max(s.UserData.maxFreq(:));
    qh=linspace(-0.5,.5,50);qk=qh; ql=qh; w=linspace(0.01,s.UserData.maxFreq*1.2,51);
    f=iData(s,[],qh,qk,ql',w);
  end
  
  if (~isfield(s.UserData,'DOS') || isempty(s.UserData.DOS)) ...
    && isfield(s.UserData,'FREQ') && ~isempty(s.UserData.FREQ)
    nmodes = size(s.UserData.FREQ,2);
    if isempty(n)
      n = max(nmodes*10, 100);
    end
    index= find(imag(s.UserData.FREQ) == 0);
    dos_e = s.UserData.FREQ(index);
    omega_e = linspace(0,max(dos_e(:))*1.2, n);
    [dos_e,omega_e]=hist(dos_e,omega_e);
    dos_factor = size(s.UserData.FREQ,2) / trapz(omega_e(:), dos_e(:));
    dos_e = dos_e * dos_factor ; % 3n modes per unit cell
    DOS=iData(omega_e,dos_e);
    DOS.Title = [ 'Total DOS ' s.Name ];
    DOS.Label = '';
    xlabel(DOS,'Energy [meV]'); 
    ylabel(DOS,[ 'Total DOS/unit cell ' strtok(s.Name) ]);
    DOS.Error=0; s.UserData.DOS=DOS;
    % partial phonon DOS (per mode) when possible
    pDOS = [];
    for mode=1:nmodes
      f1 = s.UserData.FREQ(:,mode);
      index = find(imag(f1) == 0);
      dos_e = hist(f1(index),omega_e, n);
      dos_e = dos_e * dos_factor ; % normalize to the total DOS
      DOS=iData(omega_e,dos_e);
      DOS.Title = [ 'Mode [' num2str(mode) '] DOS ' s.Name ]; 
      DOS.Label = '';
      xlabel(DOS,'Energy [meV]'); 
      ylabel(DOS,[ 'Partial DOS[' num2str(mode) ']/unit cell ' strtok(s.Name) ]);
      DOS.Error = 0;
      pDOS = [ pDOS DOS ];
    end
    s.UserData.DOS_partials=pDOS;
    clear f1 index dos_e omega_e dos_factor DOS pDOS
  end
  
  % compute the gDOS from the the S(q,w) integral. A pre-factor is missing ?
  if ~isempty(f) && isa(f, 'iData') && ndims(f) == 4 ...
    && (~isfield(s.UserData,'gDOS') || isempty(s.UserData.gDOS))
    gDOS = camproj(f, 4);
    s.UserData.gDOS = [];
    s.UserData.gDOS = gDOS;
  end
  
  if ~isempty(inputname(1))
    assignin('caller',inputname(1),s);
  end

  % get the DOS and other output
  if isfield(s.UserData,'DOS') && ~isempty(s.UserData.DOS)
    DOS = s.UserData.DOS;
  end
  if isfield(s.UserData,'DOS_partials') && numel(s.UserData.DOS_partials) > 0
    DOS_partials = s.UserData.DOS_partials;
  end
  if isfield(s.UserData,'gDOS') && ~isempty(s.UserData.gDOS)
    gDOS = s.UserData.gDOS;
    if isempty(DOS)
      DOS = gDOS;
      gDOS=[];
    end  
  end
   
