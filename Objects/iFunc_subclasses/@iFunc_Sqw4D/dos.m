function [DOS, DOS_partials] = dos(s, n)
% iFunc_Sqw4D: dos: compute the density of states (vDOS)
%
%  The routine can be used for 4D models.
%    when used on 4D models S(HKL,w), the vDOS is computed.
%
%    DOS = dos(s)    returns the vibrational density of states (vDOS)
%      the vDOS and the partials per mode are also stored in the UserData.
%    DOS = dos(s, n) does the same with n-bins on the vDOS (n=100)
%
%    When the DOS has already been computed, it is used as is. To force a
%    recomputation, specify a different number of bins 'n' or set:
%      s.UserData.DOS=[];
%    or use DOS = dos(s,'force')
%
%    To smooth the resulting distribution, use:
%      sDOS = smooth(DOS); plot(sDOS);
%
% input:
%   s: S(q,w) 4D model (iFunc_Sqw4D)
%   n: number of energy values (integer). Optional. Default is to nmodes*10
%
% output:
%   DOS:   DOS(w)   (1D iData versus energy)
%
% Example: Sqw=sqw_cubic_monoatomic; D=dos(Sqw);
% (c) E.Farhi, ILL. License: EUPL.

% conventions:
% omega = Ei-Ef = energy lost by the neutron
%    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%    omega < 0, neutron gains energy, anti-Stokes

  DOS=[]; DOS_partials=[];
  if nargin == 0, return; end
  if nargin < 2, n = []; end
  
  % handle array of objects
  if numel(s) > 1
    for index=1:numel(s)
      DOS = [ DOS feval(mfilename, s(index), n) ];
    end
    if ~isempty(inputname(1))
      assignin('caller',inputname(1),s);
    end
    return
  end

  % compute
  [DOS, DOS_partials, s] = sqw_phonon_dos_4D(s, n);
  if ~isempty(inputname(1))
    assignin('caller',inputname(1),s);
  end
  
  % plot
  if nargout == 0 && ~isempty(DOS)
    fig=figure; 
    DOS = s.UserData.DOS;
    xlabel(DOS,[ 'Energy' ]);
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
    set(fig, 'NextPlot','new');
  end
  
% ------------------------------------------------------------------------------

function [DOS, DOS_partials, s] = sqw_phonon_dos_4D(s, n)
  % sqw_phonon_dos_4D: compute the phonon/vibrational density of states (vDOS)
  %
  % input:
  %   s: phonon S(q,w) (4D) Model or Data set (iFunc/iData)
  %
  % output:
  %   DOS:          phonon/vibrational density of states (iData)
  %   DOS_partials: phonon/vibrational density of state partials (iData array)
  %
  % Example: DOS=sqw_phonon_dos(sqw_cubic_monoatomic('defaults'))
  % (c) E.Farhi, ILL. License: EUPL.
  
  % NOTE: for partial DOS per atom 'i', see Reichardt Eq (2.9) p 4
  %   see also: Schober JNR 2014 Eq (9.105)
  %   g_i(hw) = sum_j |e_i,j|^2 delta(hw - hw_j(q))
  % i.e. we weight the DOS with polarisation vector norms per atom

  DOS = []; f=[]; DOS_partials = [];
  
  % must be 4D iFunc or iData
  if ~nargin
    return
  end
  
  if nargin < 2, n=[]; end
  if strcmp(n, 'force')
    s.UserData.DOS = [];
  end
  
  % first get a quick estimate of the max frequency
  if  ~isfield(s.UserData,'DOS') || isempty(s.UserData.DOS) || (~isempty(n) && prod(size(s.UserData.DOS)) ~= n)
    maxFreq = max(s);
    
    % evaluate the 4D model onto a mesh filling the Brillouin zone [-0.5:0.5 ]
    s.UserData.DOS     = [];  % make sure we re-evaluate again on a finer grid
    qh=linspace(-0.5,.5,50);qk=qh; ql=qh; w=linspace(0.01,maxFreq*1.2,51);
    f=iData(s,[],qh,qk,ql',w);
    % force to evaluate on a finer grid
    if ~isfield(s.UserData,'FREQ') || isempty(s.UserData.FREQ)
      qk=linspace(0,0.5,30); qh=qk; ql=qk; 
      w =linspace(0.01,maxFreq*1.2,11);
      f =iData(s,[],qh,qk,ql',w);
    end
  end
  
  if (~isfield(s.UserData,'DOS') || isempty(s.UserData.DOS)) ...
    && isfield(s.UserData,'FREQ') && ~isempty(s.UserData.FREQ)
    nmodes = size(s.UserData.FREQ,2);
    if isempty(n)
      n = max(nmodes*10, 100);
    end
    index= find(imag(s.UserData.FREQ) == 0);
    dos_e = s.UserData.FREQ(index);
    omega_e = linspace(min(dos_e(:)),max(dos_e(:))*1.2, n);
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
  elseif ~isfield(s.UserData,'FREQ') || isempty(s.UserData.FREQ)
    error([ mfilename ': Can not compute the density of states as the bare frequencies are not available (UserData.FREQ)' ]);
  end
  
  if ~isempty(inputname(1))
    assignin('caller',inputname(1),s);
  end

  % get the DOS and other output
  if isfield(s.UserData,'DOS') && ~isempty(s.UserData.DOS)
    DOS = s.UserData.DOS;
    if isfield(s.UserData,'properties')
      DOS.UserData.properties = s.UserData.properties;
    end
  end
  if isfield(s.UserData,'DOS_partials') && numel(s.UserData.DOS_partials) > 0
    DOS_partials = s.UserData.DOS_partials;
  end
   
