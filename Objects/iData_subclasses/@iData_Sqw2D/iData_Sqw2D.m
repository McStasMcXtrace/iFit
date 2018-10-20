classdef iData_Sqw2D < iData
  % iData_Sqw2D: create a 2D S(q,w) data set (iData flavour)
  %
  % The iData_Sqw2D class is a 2D data set holding a S(q,w) dynamic structure factor
  %   aka scattering function/law.
  %   The first  axis (rows)    is the Energy   transfer [meV].
  %   The second axis (columns) is the Momentum transfer [Angs-1] (wavevector). 
  %
  % This quantity usually derives from the double differential neutron cross section
  %
  %    d2 sigma 
  %   ---------   = N sigma /4pi kf/ki S(q,w)
  %   dOMEGA dEf
  %
  % conventions:
  % w = omega = Ei-Ef = energy lost by the neutron [meV]
  %    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
  %    omega < 0, neutron gains energy, anti-Stokes
  %
  % Example: s=iData_Sqw2D('SQW_coh_lGe.nc')
  %
  % Useful methods for this iData flavour:
  %
  % methods(iData_Sqw2D)
  %   all iData methods can be used.
  % iData_Sqw2D(s)
  %   convert input [e.g. a 2D iData object] into an iData_Sqw2D to give access to
  %   the methods below.
  % Sqw = ddcs2Sqw(s)
  %   convert a double differential neutron scattering cross section to a S(q,w) [Kf/Ki]
  % ddcs = Sqw2ddcs(s)
  %   convert a S(q,w) to a double differential neutron scattering cross section [Kf/Ki]
  % spw = qw2phiw(s, lambda)
  %   Convert a S(q,w) into a S(phi,w) iData (scattering angle)
  % sqt = qw2qt(s, lambda)
  %   Compute S(q,tof) from S(q,w) for given wavelength [Angs]
  % [spt,spc] = qw2phi(s, lambda)
  %   Compute S(phi,tof) from S(q,w) for given wavelength [Angs]. 
  %   Also return the S(phi,channel) as 2nd arg.
  % sab = Sab(s, M, T)
  %   Compute S(alpha,beta) from S(q,w) for given mass and temperature
  %
  % d   = dos(s)
  %   Compute the vibrational density of states.
  %
  % t   = thermochemistry(s)
  %   Compute and display thermochemistry quantities from the density of states.
  %
  % m   = moments(s)
  %   Compute the S(q,w) moments/sum rules (harmonic frequencies).
  %
  % sym = symmetrize(s)
  %   Extend the S(|q|,w) in both energy sides.
  %
  % sb  = Bosify(s)
  %   Apply the 'Bose' factor (detailed balance) to a classical data set.
  %
  % s   = deBosify(sb)
  %   Remove Bose factor (detailed balance) from an 'experimental/quantum' data set.
  %
  % p   = parseparams(s)
  %   Search for physical quantities in S(q,w) data set.
  %
  % sab = Sab(s)
  %   Convert to an S(alpha,beta) suitable for nuclear data bases (ENDF, etc).
  %
  % xs  = scattering_cross_section(s)
  %   Compute the total integrated scattering cross section
  %
  % dr  = dynamic_range(s, Ei, angles)
  %   Compute the dynamic range restricted to given incident energy and detector angles
  %
  % sq  = structure_factor(s)
  %   Compute the structure factor S(q)
  %
  % [inc, multi] = incoherent(s, q, T, m, n)
  %   Compute an estimate of the incoherent neutron scattering law in the gaussian approximation (Sjolander)
  %
  % [coh] = coherent(inc, sq)
  %   Compute an estimate of the coherent S(q,w) from an incoherent S(q,w) and a structure factor (Skold)
  %
  % [gDOS,M]     = multi_phonons(ss)
  %   Compute the integrated multi-phonon DOS terms from an initial density of states (Sjolander)
  %
  % saveas(s, filename, 'McStas'|'Sqw'|'inx'|'spe')
  %   Save the S(q,w) as a McStas Sqw, INX or ISIS SPE file format
  %
  % input:
  %   can be a 2D iData or filename to generate a 2D Sqw object.
  %
  % output: an iData_Sqw2D object
  %
  % See also: iData, iData_Sab, iData_vDOS, iFunc_Sqw2D
  % (c) E.Farhi, ILL. License: EUPL.

  properties
  end
  
  methods
    % main instantiation method
    function obj = iData_Sqw2D(s)
      % iData_Sqw2D: create the iData_Sqw2D subclass
      %
      %   convert: anything -> iData_Sqw2D
      %
      % input:
      %   must be a 2D iData object holding a S(q,w), S(phi,t), S(alpha,beta) or S(phi,w)
      %   or file name.
      %
      % Example: s=iData_Sqw2D('SQW_coh_lGe.nc');
      obj = obj@iData;
      obj.class = mfilename;
      
      if ~nargin, return; end  % empty object
      
      % convert/test
      if     isa(s, mfilename)   m = s;
      elseif isa(s, 'iData_Sab') m = Sab2Sqw(s);
      else
        m = Sqw_check(s, 'qw'); % check and possibly convert to Sqw
        if ~isa(m, 'iData') || any(isempty(m)) || any(ndims(m) ~= 2)
          error([ mfilename ': the given input ' class(s) ' does not seem to be convertible to iData_Sqw2D.' ])
        end
      end
      
      % copy all properties
      obj = copy_prop(obj, m);
      obj = commandhistory(obj, mfilename, s);
      label(obj, 0, [  mfilename '(' label(obj, 0) ')' ]);
 
    end % iData_Sqw2D constructor
    
    % parameters (search for parameters in iData)
    function parameters = parseparams(s)
      % iData_Sqw2D: parseparams: search for physical quantities in object.
      % This search is also done when creating iData_Sqw2D objects.
      %
      %   iData_Sqw2D -> physical parameters
      [s,parameters,fields] = Sqw_parameters(s);
      if length(inputname(1))
        assignin('caller',inputname(1),s);
      end
    end
    
    function f = iData(self)
      % iData_Sqw2D: iData: convert a iData_Sqw2D back to iData
      %
      % convert: iData_Sqw2D -> iData
      f = [];
      for index=1:numel(self)
        f1   = copy_prop(iData, self(index));
        f1   = commandhistory(f1, 'iData', self(index));
        label(f1, 0, [  'iData' '(' label(self(index), 0) ')' ]);
        f = [ f f1 ];
      end
    end
    
    % sound_velocity ?
    % MSD ?
    % diffusion constant ?
    % compressibility ?
    % g(r) pdf
    
    function [inc, single, multi] = incoherent(s, varargin)
      % iData_Sqw2D: incoherent: incoherent neutron scattering law estimate in the incoherent gaussian approximation
      %
      %   [inc, single, multiple] = incoherent(s, q, T, m, n, DW)
      %
      % compute: iData_Sqw2D -> generalised Density of States -> Incoherent approximation S(q,w)
      %
      % The result is the dynamic structure factor (scattering law) for neutrons, in
      %   the incoherent gaussian approximation. If you wish to obtain the intermediate
      %   scattering function I(q,t), use [inc, Iqt] = incoherent(dos(s), ...)
      %   These should be e.g. multiplied by the neutron scattering bound cross 
      %   section 'sigma_inc' [barns]. This calculation includes the Debye-Waller factor.
      %
      % This implementation is in principle exact for an isotropic monoatomic material,
      %   e.g. a liquid or powder.
      % This methodology is equivalent to the LEAPR module of NJOY ("phonon expansion")
      %   to compute S(alpha,beta) from a vibrational density of states.
      %
      % For a poly-atomic material with a set of non-equivalent atoms with relative 
      %   concentration Ci, mass Mi and bound scattering cross section sigma_i, 
      %   one should use:
      %
      %   sigma = sum_i Ci sigma_i                              weighted cross section
      %   m     = [sum_i Ci sigma_i]/[sum_i Ci sigma_i/Mi]      weighted mass
      %
      % conventions:
      % w = Ei-Ef = energy lost by the neutron
      %    w > 0, neutron looses energy, can not be higher than Ei (Stokes)
      %    w < 0, neutron gains energy, anti-Stokes
      %
      % Reference:
      %   H. Schober, Journal of Neutron Research 17 (2014) 109–357
      %     DOI 10.3233/JNR-140016 (see esp. pages 328-331)
      %   V.S. Oskotskii, Sov. Phys. Solid State 9 (1967), 420.
      %   A. Sjolander, Arkiv for Fysik 14 (1958), 315.
      %
      % syntax:
      %   [inc, single, multi] = incoherent(sqw)
      %   [inc, single, multi] = incoherent(sqw, q, T, m, n, DW)
      %   [inc, single, multi] = incoherent(sqw, 'q', q, 'T', T, 'm', m, 'n', n, 'DW', dw)
      %
      % Missing arguments (or given as [] empty), are searched within the initial 
      %   data set. Input arguments can be given in order, or with name-value 
      %   pairs, or as a structure with named fields.
      %
      % input:
      %   sqw:the scattering law S(q,w) [iData_Sqw2D]
      %   q:  the momentum axis [Angs-1, vector]
      %   T:  temperature [K]
      %   m:  mass of the scattering unit [g/mol]
      %   n:  number of iterations in the series expansion, e.g. 5
      %   DW: Debye-Waller coefficient gamma=<u^2> [Angs^2] e.g. 0.005
      %       The Debye-Waller function is      2W(q)=gamma*q^2
      %       The Debye-Waller factor   is exp(-2W(q))
      %
      % output:
      %   inc:    total S(q,w), single+multi-phonons [iData_Sqw2D]
      %   single: so-called single-phonon incoherent S(q,w) [iData_Sqw2D]
      %   multi:  so-called multi-phonon incoherent S(q,w)  [iData_Sqw2D]
      %
      % Example:
      %   s   = iData_Sqw2D('D2O_liq_290_coh.sqw.zip');
      %   inc = incoherent(s); plot(log10(inc));
      %
      % See also: iData_Sqw2D/multi_phonons_dos
      % (c) E.Farhi, ILL. License: EUPL.
      g   = dos(s, varargin{:});
      inc = incoherent(g, varargin{:});
      multi = plus(inc(3:end)); % multi-phonon
      single= plus(inc(1:2));   % elastic+1-phonon
      inc   = plus(inc);        % total
      % update Command
      inc   = commandhistory(inc,    'incoherent', s, varargin{:});
      multi = commandhistory(multi,  'incoherent', s, varargin{:});
      single= commandhistory(single, 'incoherent', s, varargin{:});
      label(inc,    0, [  'incoherent' '(' label(s, 0) ')' ]);
      label(multi,  0, [  'multi-phonons' '(' label(s, 0) ')' ]);
      label(single, 0, [  'single-phonon' '(' label(s, 0) ')' ]);
      if nargout == 0
        fig=figure; 
        h  =subplot(log10([inc single multi]),'view2'); 
        set(fig, 'NextPlot','new');
      end
    end % incoherent
    
    function [coh] = coherent(self, sq)
      % iData_Sqw2D: coherent: compute the coherent scattering cross section from the incoherent and structure factor, in the Skold approximation
      %
      %   coh = coherent(inc, sq)
      %
      % compute: iData_Sqw2D incoherent S(q,w) + S(q) -> coherent S(q,w)
      %
      % The Skold approximation is:
      %
      %   Scoh(q,w) = Sinc(q/sqrt(S(q)), w) S(q)
      %
      % The result should be e.g. multiplied by the neutron scattering bound cross 
      % section 'sigma_coh' [barns].
      %
      % Reference: K. Skold, Phys. Rev. Lett. 19, 1023 (1967).
      %
      % input:
      %   inc: incoherent S(q,w) [iData_Sqw2D]
      %   sq:  S(q)              [double or iData]
      %        can also be a single interatomic distance d-spacing
      % output:
      %   coh: coherent estimate [iData_Sqw2D]
      
      if nargin < 2, sq = []; end
      
      % make sure we use a regular grid
      inc = meshgrid(self,'vector');
      q   = getaxis(inc, 2);
      
      if isempty(sq)
        sq = [ 3 .2]; % will use Percus-Yevick
      end
      if prod(size(sq)) <= 2
        % assume sq is a mean interatomic distance -> Percus-Yevick model
        if isscalar(sq), p = [ double(sq) .2 ];
        else p =double(sq); end
        sq = sf_hard_spheres(p);
      end
      if isa(sq, 'iFunc')
        sq = feval(sq, q);
      end
      if isa(sq, 'iData') && ndims(sq)==1
        sq = interp(sq, q); % make sure we use same q values in S(q,w) and S(q)
      end
      
      sq = double(sq); sq = sq(:);
      
      % K. Skold, Phys. Rev. Lett. 19, 1023 (1967).
      % Scoh(q,w) = Sinc(q/sqrt(S(q)), w) S(q)
      
      % index method: we directly use index values to avoid interpolations
      % index of new Skold q value in initial Sinc data set
      
      % number of q values in initial and final data sets
      nq = size(self, 2); 
      index_q_skold = round( (1:nq)./ sqrt(sq') );
      index_q_skold(index_q_skold < 1)  = 1;
      index_q_skold(index_q_skold > nq) = nq;
      index_q_skold(~isfinite(index_q_skold)) = [];
      
      coh = copyobj(inc);
      signal_inc = getaxis(inc, 0);
      signal_coh = 0*signal_inc;
      signal_coh(:,1:numel(index_q_skold)) = signal_inc(:,index_q_skold);  % Sinc(q/sqrt(sq),w)
      coh = set(coh, 'Signal', signal_coh);
      
      % interpolation method: compute the new Sinc with modified q axis
      % qinc = getaxis(inc, 2)./sqrt(sq);
      % sinc = interp(s, qinc, w).*sq;
      % interpolate on the initial grid
      % coh  = interp(coh, q,w);
      
      % get the coherent estimate
      coh  = coh.*sq';

    end % coherent
    
    function [G,multi,g] = multi_phonons(s, varargin)
      % iData_Sqw2D: multi_phonons: compute the integrated multi-phonon intensity from a scattering law
      %
      % compute: iData_Sqw2D -> generalised Density of States -> multi-phonons DOS
      %
      % output:
      %   G:      total gDOS [sum p=1..Inf]
      %   multi:  multi-phonon gDOS contribution [p=2...]
      %   g:      gDOS [p=1]
      g     = dos(s);
      G     = multi_phonons(g, varargin{:});
      multi = plus(G(2:end));
      G     = plus(G);
      G     = commandhistory(G,    'multi_phonons', s, varargin{:});
      multi = commandhistory(multi,'multi_phonons', s, varargin{:});
      label(G,     0, [  'total gDOS' '(' label(s, 0) ')' ]);
      label(multi, 0, [  'multi_phonons gDOS' '(' label(s, 0) ')' ]);
      if nargout == 0
        fig=figure; 
        h  =plot([ G multi ]); 
        set(fig, 'NextPlot','new');
      end
    end
    
    function g = gdos(self, varargin)
      % iData_Sqw2D: gdos: compute the generalised density of states (gDOS) from a S(q,w)
      %
      % compute: iData_Sqw2D -> generalised Density of States
      %
      % See: iData_Sqw2D/dos
      g = dos(self, varargin{:});
    end
    
    function v = vdos(self, varargin)
      % iData_Sqw2D: vdos: compute the 'true' vibrational density of states
      %   by removing the multi-phonons from the initial estimate from the 
      %   Bredov/Oskotskii formula.
      %
      % compute: iData_Sqw2D -> generalised Density of States -> vibrational DOS
      
      % we first compute an initial estimate of the vDOS
      g = dos(self);
      v = vdos(g, varargin{:});  % then get the true vDOS using the iData_vDOS method.
    end
    
    function spw = qw2phiw(self, varargin)
      % iData_Sqw2D: qw2phiw: convert a S(q,w) into a S(phi,w) iData (scattering angle)
      %
      % convert: iData_Sqw2D S(q,w) -> S(phi,w)
      %
      % spw = qw2phiw(self)
      % spw = qw2phiw(self, lambda)
      
      spw = Sqw_q2phi(self, varargin{:});
      spw = commandhistory(spw, 'qw2phiw', self, varargin{:});
      spw.Label = 'S(phi, w)';
      label(spw, 0, [  'qw2phiw' '(' label(self, 0) ')' ]);
      if nargout == 0
        fig=figure; 
        h  =plot(log10(spw)); 
        set(fig, 'NextPlot','new');
      end
    end
    
    function sqt = qw2qt(self, varargin)
      % iData_Sqw2D: qw2qt: convert a S(q,w) into a S(q,tof) iData (time-of-flight from sample)
      %
      % convert: iData_Sqw2D S(q,w) -> S(q,tof)
      %
      % sqt = qw2qt(self)
      % sqt = qw2qt(self, lambda)
      sqt = Sqw_e2t(self, varargin{:});
      sqt = commandhistory(sqt, 'qw2qt', self, varargin{:});
      sqt.Label = 'S(q, tof)';
      label(sqt, 0, [  'qw2qt [e2t]' '(' label(self, 0) ')' ]);
      if nargout == 0
        fig=figure; 
        h  =plot(log10(sqt)); 
        set(fig, 'NextPlot','new');
      end
    end
    
    function [spt, spc, spw] = qw2phi(self, varargin)
      % iData_Sqw2D: qw2phi: convert a S(q,w) into S(phi,tof) iData (scattering angle,ToF)
      %   This method returns S(phi,tof), S(phi,tof channels) and S(phi,w)
      %
      % convert: iData_Sqw2D S(q,w) -> [S(phi,tof) S(phi,tof channels) S(phi,w)]
      %
      %   [spt, spc, spw] = qw2phi(sqw);
      %   [spt, spc, spw] = qw2phi(sqw, lambda);
      %
      % input:
      %   sqw:    S(q,w) as an iData_Sqw2D object
      %   lambda: wavelength used for conversion [Angs]. When not given , it is
      %           searched in the S(q,w) data set.
      %
      % output:
      %   spt:    S(phi,tof)      with phi=scattering angle [deg], t=time-of-flight from sample
      %   spc:    S(phi,channel)  with the time-of-flight axis given in time channels
      %   spw:    S(phi,w)        with phi=scattering angle [deg]
      spw        = Sqw_q2phi(self, varargin{:});
      [spt, spc] = Sqw_e2t(spw, varargin{:});
      if nargout == 0
        fig=figure; 
        h  =plot(log10(spt)); 
        set(fig, 'NextPlot','new');
      end
    end
    
    function ddcs = Sqw2ddcs(s, lambda, inverse)
      % iData_Sqw2D: Sqw2ddcs: convert a S(q,w) into a double differential cross-section (q,w)
      %   i.e. multiply by Kf/Ki
      %
      % convert: iData_Sqw2D S(q,w) -> d2(sigma)/dOmega/dE = N.sigma Kf/Ki S(q,w)
      %
      %   s = Sqw2ddcs(s, lambda)
      %
      % input:
      %   s:      S(q,w) as an iData_Sqw2D object
      %   lambda: wavelength used for conversion [Angs]. When not given , it is
      %           searched in the S(q,w) data set.
      %
      % output:
      %   s:      double differential cross-section [Kf/Ki S(q,w)]
      
      ddcs = [];
      if isempty(s), return; end
      if nargin < 2, lambda = []; end
      if nargin < 3, inverse = false; else inverse = true; end
      if isempty(lambda)
        [s,lambda,distance,chwidth] = Sqw_search_lambda(s);
      else
        [s,~,distance,chwidth] = Sqw_search_lambda(s);
      end
      if isempty(lambda)
        lambda   = 2.36;
        disp([ mfilename ': ' s.Tag ' ' s.Title ' using <wavelength>               =' num2str(lambda) ' [Angs]']);
      end

      SE2V = 437.393377;        % Convert sqrt(E)[meV] to v[m/s]
      V2K  = 1.58825361e-3;     % Convert v[m/s] to k[1/AA]
      K2V  = 1/V2K;
      VS2E = 5.22703725e-6;     % Convert (v[m/s])**2 to E[meV]
      
      Ki   = 2*pi/lambda;
      Vi   = K2V*Ki;
      Ei   = VS2E*Vi.^2;
      % compute final energy
      hw   = getaxis(s, 1);
      Ef   = Ei - hw;
      kikf = sqrt(Ei./Ef);
      if inverse
        kikf = 1./kikf; % DDCS -> Sqw
      end
      
      if ~inverse
        % check if initial data set is classical: apply Bose factor
        if (isfield(s,'classical') || ~isempty(findfield(s, 'classical')))
          classical = get(s,'classical');
        else classical = false;
        end
        if ~isempty(classical) && classical(1)
          s = Bosify(s);
        end
      end
      
      ddcs = copyobj(s) .* kikf;
      setalias(ddcs, 'IncidentWavelength', lambda);
      
      if inverse
        ddcs = commandhistory(ddcs, 'ddcs2Sqw', s, lambda, inverse);
        ddcs.Label = 'S(q, w)';
        label(ddcs, 0, [  'ddcs2Sqw' '(' label(s, 0) ')' ]);
      else
        ddcs = commandhistory(ddcs, 'Sqw2ddcs', s, lambda, inverse);
        ddcs.Label = 'DDCS(q, w)';
        label(ddcs, 0, [  'DDCS' '(' label(s, 0) ')' ]);
      end
    
    end % Sqw2ddcs
    
    function self = ddcs2Sqw(ddcs, varargin)
      % iData_Sqw2D: ddcs2Sqw: convert a double differential cross-section (q,w) into a S(q,w)
      %   i.e. divide by Kf/Ki
      %
      % convert: iData_Sqw2D d2(sigma)/dOmega/dE -> S(q,w)
      %
      %   s = ddcs2Sqw(s, lambda)
      %
      % input:
      %   s:      double differential cross-section [Kf/Ki S(q,w)]
      %   lambda: wavelength used for conversion [Angs]. When not given , it is
      %           searched in the S(q,w) data set.
      %
      % output:
      %   s:      S(q,w) as an iData_Sqw2D object
      %
      %   s = ddcs2Sqw(s, lambda)
      %
      %   d2(sigma)/dOmega/dE = N.sigma Kf/Ki S(q,w)
      
      if nargin > 1
        self = Sqw2ddcs(ddcs, varargin{:}, 'inverse');
      else
        self = Sqw2ddcs(ddcs, [], 'inverse');
      end
    
    end % qw2ddcs
    
    function s = Sab(self, varargin)
      %  iData_Sqw2D: Sab: convert a 2D S(q,w) into an S(alpha,beta). 
      %
      % convert: iData_Sqw2D S(q,w) -> S(alpha,beta)
      %
      %  syntax: sab = Sab(sqw, M, T)
      %
      %  The S(alpha,beta) is a representation of the dynamic structure factor 
      %  using unitless momentum and energy variables defined as:
      %     alpha= h2q2/2MkT  = (Ei+Ef-2*mu*sqrt(Ei*Ef))/AkT
      %     beta = -hw/kT     = (Ef-Ei)/kT
      %     A    = M/m
      %     mu   = cos(theta) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf)
      %  
      % input:
      %   s:  S(q,w) data set e.g. 2D data set with q [Angs-1] as 1st axis (rows), energy [meV] as 2nd axis.
      %   M:  molar weight of the atom/molecule in [g/mol].
      %     when omitted or empty, it is searched as 'weight' or 'mass' in the object.
      %   T: when given, Temperature to use. When not given or empty, the Temperature
      %      is searched in the object. The temperature is in [K]. 1 meV=11.605 K.
      % output:
      %   sab: S(alpha,beta) 2D data set, classical (iData_Sab)
      %
      % conventions:
      % w = omega = Ei-Ef = energy lost by the neutron [meV]
      %    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
      %    omega < 0, neutron gains energy, anti-Stokes
      %
      % references: M. Mattes and J. Keinert, IAEA INDC(NDS)-0470, 2005.
      %             R. E. MacFarlane, LA-12639-MS (ENDF-356), 1994.
      %
      % Example: sqw = iData_Sqw2D('SQW_coh_lGe.nc');
      %          sab = Sab(sqw,72.6,1235);
      %          plot(log(sab)
      %
      % (c) E.Farhi, ILL. License: EUPL.
      s = Sqw2Sab(self, varargin{:});  % private
    end
    
    function f = saveas(self, varargin)
      % iData_Sqw2D: saveas: save S(q,w) into a file.
      %
      % convert: iData_Sqw2D S(q,w) -> file
      %
      % syntax: saveas(sqw2D, filename, format)
      %
      %   with format as mcstas, sqw, spe, inx or any other iData supported format.
      %
      % all iData.saveas formats are available, and in addition:
      %
      %   McStas (for Isotropic_Sqw component)
      %   INX (INX/MuPhoCor)
      %   SPE (MSlice)
      %
      if numel(varargin) >= 2
        switch lower(varargin{2})
        case {'mcstas','sqw'}
          f = Sqw_McStas(self, varargin{1});
        case 'inx'
          spe = Sqw_q2phi(self);  % S(phi, w)
          f = write_inx(spe, varargin{1});
        case 'spe'
          spe = Sqw_q2phi(self);  % S(phi, w)
          f = write_spe(spe, varargin{1});
        otherwise
          f = saveas@iData(self, varargin{:});
        end
      else
        f = saveas@iData(self, varargin{:});
      end
      
    end
    
  end
  
end
