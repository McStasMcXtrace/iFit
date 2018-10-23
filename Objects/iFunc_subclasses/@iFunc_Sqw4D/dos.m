function [DOS, DOS_partials] = dos(s, n, nQ, method)
% iFunc_Sqw4D: dos: compute the density of states (vDOS)
%
%  The routine can be used with 4D models to compute the vibrational density of 
%   states (vDOS), aka phonon spectrum.
%
%    DOS = dos(s)    returns the vibrational density of states (vDOS)
%      the vDOS and the partials per mode are also stored in the model UserData.
%    DOS = dos(s, n) 
%      does the same with n-bins on the vDOS (n=100). 
%      'n' can also be a vector of energy values.
%    DOS = dos(s, n, nQ) 
%      does the same with nQ-bins for the HKL average (n=21)
%    DOS = dos(s, n, nQ, method) 
%      specifies the method among '4d', 'powder' and 'fast' estimate.
%
%    To smooth the resulting distribution, use:
%      sDOS = smooth(DOS); plot(sDOS);
%
%    If the DOS has already been calculated, it is re-used. To force calculation
%    use a non-0 value or 'force' for 'n'.
%
% input:
%   s:      S(q,w) 4D model (iFunc_Sqw4D)
%   n:      number/vector of energy values (scalar/vector). Optional. Default is nmodes*10
%   nQ:     number of Q-grid binning (integer). Optional. Default is 21
%   method: optional, can be '4D' (default) or 'powder' or 'fast'
%
% output:
%   DOS:          DOS(w)   (1D iData versus energy)
%   DOS_partials: DOS(w) per mode (1D iData array, versus energy)
%
% See also: iFunc_Sqw4D/max
%
% Example: Sqw=sqw_cubic_monoatomic; D=dos(Sqw);
% (c) E.Farhi, ILL. License: EUPL.

% conventions:
% omega = Ei-Ef = energy lost by the neutron
%    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%    omega < 0, neutron gains energy, anti-Stokes

  DOS=[]; DOS_partials=[];
  if nargin == 0, return; end
  if nargin < 2,  n  = []; end
  if nargin < 3,  nQ = []; end
  if nargin < 4,  method=''; end
  
  if strcmp(n, 'force'), n=0; end
  if ischar(n)  && isempty(method), method = n;  n=[]; end
  if ischar(nQ) && isempty(method), method = nQ; nQ=[]; end
  if isempty(method), method = '4D'; end
  
  % handle array of objects
  if numel(s) > 1
    for index=1:numel(s)
      DOS = [ DOS feval(mfilename, s(index), n, nQ, method) ];
    end
    if ~isempty(inputname(1))
      assignin('caller',inputname(1),s);
    end
    return
  end
  
  if isfield(s.UserData,'DOS') && ~isempty(s.UserData.DOS)
    DOS = s.UserData.DOS;
  end
  if isfield(s.UserData,'DOS_partials') && numel(s.UserData.DOS_partials)
    DOS_partials = s.UserData.DOS_partials;
  end
  
  if ~isempty(n) || isempty(DOS)

    % compute frequencies
    switch lower(method)
    case '4d'
      FREQ = sqw_phonon_dos_4D(s, nQ);
    case 'fast'
      FREQ=[];
    otherwise
      FREQ = sqw_phonon_dos_powder(s, nQ);
    end
    
    if isempty(FREQ)
      try
        [~,DOS] = max(s);
      end
    else
      % compute DOS
      [DOS, DOS_partials] = dos_getdos_from_FREQ(FREQ, n, s.Name);
      clear FREQ
    end

    % transfer UserData stuff
    for f={'properties','calc','configuration','options','FORCES','dir','maxFreq'}
      if isfield(s.UserData, f{1})
        DOS.UserData.(f{1}) = s.UserData.(f{1});
      end
    end
    DOS = iData_vDOS(DOS);
    parameters = get(DOS, 'parameters');
    for index=1:numel(s.ParameterValues)
      if isfinite(s.ParameterValues(index))
        parameters.(strtok(s.Parameters{index})) = s.ParameterValues(index);
      end
    end
    setalias(DOS, 'parameters', parameters);
    DOS.UserData.ModelParameters = parameters;
    
    if ~isempty(inputname(1))
      s.UserData.DOS          = DOS;
      s.UserData.DOS_partials = DOS_partials;
      assignin('caller',inputname(1),s);
    end
  end % else re-use existing stored DOS
  
  % plot
  if nargout == 0 && ~isempty(DOS)
    fig=figure; 
    xlabel(DOS,[ 'Energy' ]);
    % plot any partials first
    h=plot(DOS_partials);
    if iscell(h), h=[ h{:} ]; end
    set(h,'LineStyle','--');
    hold on

    % plot total DOS and rotate
    h=plot(DOS); set(h,'LineWidth',2);
    set(fig, 'NextPlot','new');
  end

% ------------------------------------------------------------------------------
function FREQ = sqw_phonon_dos_4D(s, nQ)
  % sqw_phonon_dos_4D: compute the phonon/vibrational density of states (vDOS)
  %
  % input:
  %   s: phonon S(q,w) (4D) Model or Data set (iFunc/iData)
  %
  % output:
  %   FREQ:         mode frequencies
  %
  % Example: DOS=sqw_phonon_dos(sqw_cubic_monoatomic('defaults'))
  % (c) E.Farhi, ILL. License: EUPL.

  FREQ = [];
  if isempty(nQ) || nQ <= 0, nQ=21; end
  
  qh = linspace(-0.5,0.5,nQ);qk=qh; ql=qh; 
  w  = linspace(0.01,100,5);
  try
    f=feval(s,[],qh,qk,ql',w);  % this updates the object with mode frequencies
    clear f
  end
  
  % get bare frequencies
  if isfield(s.UserData, 'FREQ')
    FREQ = s.UserData.FREQ;
  end
  
% ------------------------------------------------------------------------------
function FREQ = sqw_phonon_dos_powder(s, nQ)
  
  FREQ = [];
  if isempty(nQ) || nQ <= 0, nQ=41; end
  
  % create a 4D hklw space grid
  qx=[]; qy=qx; qz=qx;
  w  = linspace(0.01,100,5); % no need for this energy axis, we use bare FREQ
  n=100; % n points at constant |Q|
  x = linspace(-0.5, 0.5, nQ);  % in [rlu]
  for index=1:numel(x)
    XYZ = randn(3,n);
    XYZ = x(index)*bsxfun(@rdivide,XYZ,sqrt(sum(XYZ.^2,1)));
    qx = [ qx XYZ(1,:) ];
    qy = [ qy XYZ(2,:) ];
    qz = [ qz XYZ(3,:) ];
  end
  clear XYZ

  try
    f = feval(s,[],qx,qy,qz,w);
    clear f
  end
  
  if isfield(s.UserData, 'FREQ')
    FREQ = s.UserData.FREQ;
  end
  
% ------------------------------------------------------------------------------
function [DOS, pDOS] = dos_getdos_from_FREQ(FREQ, n, Name)

  % NOTE: for partial DOS per atom 'i', see Reichardt Eq (2.9) p 4
  %   see also: Schober JNR 2014 Eq (9.105)
  %   g_i(hw) = sum_j |e_i,j|^2 delta(hw - hw_j(q))
  % i.e. we weight the DOS with polarisation vector norms per atom

  nmodes = size(FREQ,2);
  if ischar(n), n=[]; end
  if isempty(n) || n <= 0
    n = nmodes*10;
  end
  
  % compute the DOS histogram
  index           = find(imag(FREQ) == 0);
  dos_e           = FREQ(index);
  if isscalar(n)
    omega_e         = linspace(min(dos_e(:)),max(dos_e(:))*1.2, n);
    dn = 0;
  else
    omega_e = [ n(1)-1 sort(n) n(end)+1 ];  % first bin accumulates below, and last, above.
    n = numel(omega_e)-2;
  end
  dos_e           = hist(dos_e,omega_e);
  N3              = size(FREQ,2); % number of modes = 3N
  dos_factor      = N3 / trapz(omega_e(:), dos_e(:));
  dos_e           = dos_e * dos_factor ; % 3n modes per unit cell
  if n == numel(omega_e)-2
    omega_e = omega_e(2:(end-1));
    dos_e   = dos_e(2:(end-1));
  end
  
  % create the object
  DOS                   = iData(omega_e,dos_e);
  DOS.Title             = [ 'Total DOS ' Name ];
  DOS.Label             = '';
  DOS.Error             = 0;
  xlabel(DOS,'Energy [meV]'); 
  ylabel(DOS,[ 'Total DOS/unit cell ' strtok(Name) ]);
  DOS= iData_vDOS(DOS);
  
  % partial phonon DOS (per mode) when possible
  pDOS   = [];
  for mode=1:nmodes
    f1 = FREQ(:,mode);
    index = find(imag(f1) == 0);
    dos_e = hist(f1(index),omega_e);
    dos_e = dos_e * dos_factor ; % normalize to the total DOS
    if n == numel(omega_e)-2
      dos_e   = dos_e(2:(end-1));
    end
    
    d       = iData(omega_e,dos_e);
    d.Title = [ 'Mode [' num2str(mode) '] DOS ' Name ]; 
    d.Label = '';
    d.Error = 0;
    xlabel(d,'Energy [meV]'); 
    ylabel(d,[ 'Partial DOS[' num2str(mode) ']/unit cell ' strtok(Name) ]);
    d = iData_vDOS(d);
    pDOS      = [ pDOS d ];
  end
  
  s.UserData.DOS_partials= pDOS;
