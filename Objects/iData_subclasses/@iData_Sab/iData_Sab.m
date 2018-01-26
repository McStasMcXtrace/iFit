classdef iData_Sab < iData
  % iData_Sab: create a 2D S(alpha,beta) data set (iData flavour)
  %
  % The iData_Sab class is a 2D data set holding a S(alpha,beta) dynamic 
  %   structure factor aka scattering function/law.
  %   The data set axes are beta as 1st axis (rows), alpha as 2nd axis (Angs-1).
  %
  % The S(alpha,beta) is a representation of the dynamic structure factor 
  % using unitless momentum and energy variables defined as:
  %     alpha= h2q2/2MkT  = (Ei+Ef-2*mu*sqrt(Ei*Ef))/AkT
  %     beta = -hw/kT     = (Ef-Ei)/kT                  energy gained by neutron
  %     A    = M/m
  %     mu   = cos(theta) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf)
  % This representation is common in nuclear data, neutron sections (e.g. ENDF MF7).
  %
  % conventions:
  % beta = (Ef-Ei)/kT = -hw/kT  energy gained by the neutron, unitless
  %    beta < 0, neutron looses energy, down-scattering (Stokes)
  %    beta > 0, neutron gains energy,  up-scattering (anti-Stokes)
  %
  % Example: sab=iData_Sab('SQW_coh_lGe.nc', 72.6,1235)
  %
  % Useful methods for this iData flavour:
  %
  % methods(iData_Sab)
  %   all iData methods can be used.
  % iData_Sab(sab)
  %   convert input [e.g. a 2D iData object] into an iData_Sab to give access to
  %   the methods below.
  %
  % d   = dos(sab)
  %   Compute the generalized vibrational density of states (gDOS).
  %
  % t   = thermochemistry(sab,T)
  %   Compute and display thermochemistry quantities from the gDOS.
  %
  % m   = moments(sab, M, T)
  %   Compute the S(alpha,beta) moments/sum rules (harmonic frequencies).
  %
  % sym = symmetrize(sab)
  %   Extend the S(alpha,beta) in both 'beta' sides. The resulting S(a,b) is classical/symmetric in energy.
  %
  % sb  = Bosify(sab,T)
  %   Apply the 'Bose' factor (detailed balance) to a classical data set.
  %
  % s   = deBosify(sab,T)
  %   Remove Bose factor (detailed balance) from an 'experimental/quantum' data set.
  %
  % p   = parseparams(sab)
  %   Search for physical quantities in S(alpha,beta) data set.
  %
  % sqw = Sqw(sab)
  %   Convert an S(alpha,beta) to an S(q,w), which is roughly independent of the temperature.
  %
  % xs  = scattering_cross_section(sab, Ei)
  %   Compute the total integrated scattering cross section
  %
  % dr  = dynamic_range(sab, Ei, angles)
  %   Compute the dynamic range restricted to given incident energy and detector angles
  %
  % sq  = structure_factor(sab)
  %   Compute the structure factor S(q)
  %
  % [inc,multi] = incoherent(sab)
  %   Compute the incoherent neutron scattering law estimate in the incoherent 
  %     gaussian approximation, using its density of states.
  %
  % [gDOS,M]    = multi_phonons(sab)
  %   Compute the integrated multi-phonon DOS terms from an initial density of states
  %
  % input:
  %   can be an iData or filename to generate a 2D Sab object.
  %
  % output: an iData_Sab object
  %
  % See also: iData, iData_Sqw2D, iData_vDOS
  % (c) E.Farhi, ILL. License: EUPL.
  
  properties
  end
  
  methods
    % main instantiation method
    function obj = iData_Sab(s, varargin)
      % iData_Sab: create the iData_Sab subclass
      %
      % syntax:
      %   blah = iData_Sab(something)
      %
      % convert 'something' into an iData_Sab object so that specific neutron 
      % S(alpha,beta) methods can be used.
      %
      % input:
      %   s: must be a 2D iData object holding a S(q,w), S(phi,t), S(alpha,beta) or S(phi,w)
      %      or file name with such data.
      %   M: molar mass [g/mol]. Searched in the input data set when not given.
      %   T: temperature [K]. Searched in the input data set when not given.
      % output:
      %   S(alpha,beta) data set
      %
      % Example: sab=iData_Sab('SQW_coh_lGe.nc', 72.6,1235);
      %
      % See also: iData_Sqw2D
      obj = obj@iData;
      obj.class = mfilename;
      
      if ~nargin, return; end  % empty object
      
      % convert/test
      if     isa(s, mfilename)   m = s; 
      elseif isa(s, 'iData_Sqw') m = Sqw2Sab(s, varargin{:});
      else
        m = Sqw_check(s, 'ab'); 
        if ~isa(m, 'iData') || any(isempty(m)) || any(ndims(m) ~= 2)
          error([ mfilename ': the given input ' class(s) ' does not seem to be convertible to iData_Sab.' ])
        end
      end
      
      % copy all properties
      obj = copy_prop(obj, m);
 
    end % iData_Sab constructor
    
    % parameters (search for parameters in iData)
    function parameters = parseparams(s)
      % iData_Sab: parseparams(sab): search for physical quantities in object.
      %   The initial object is updated with a 'parameter' property.
      %
      % syntax:
      %   p = parseparams(sab)
      %
      % input:  sab: S(alpha,beta) [iData_Sab]
      % output: physical parameters [struct]
      %
      % Example: sab=iData_Sab('SQW_coh_lGe.nc', 72.6,1235); parseparams(sab)
      %
      % See also: iData_Sqw2D/parseparams
      [s,parameters,fields] = Sqw_parameters(s);
      if nargout == 0 && length(inputname(1))
        assignin('caller',inputname(1),s);
      end
    end
    
    function s = Sqw(self)
      % sqw = Sqw(Sab, M, T)
      %  iData_Sab: Sqw: convert a 2D S(alpha,beta) into an S(q,w).
      %
      %  The S(alpha,beta) is a representation of the dynamic structure factor 
      %  using unitless momentum and energy variables defined as:
      %     alpha= h2q2/2MkT  = (Ei+Ef-2*mu*sqrt(Ei*Ef))/AkT
      %     beta = -hw/kT     = (Ef-Ei)/kT
      %     A    = M/m
      %     mu   = cos(theta) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf)
      %  
      % input:
      %   s:  S(alpha,beta) data set e.g. 2D data set with beta as 1st axis (rows), alpha as 2nd axis (columns).
      %   M:  molar weight of the atom/molecule in [g/mol].
      %     when omitted or empty, it is searched as 'weight' or 'mass' is the object.
      %   T: when given, Temperature to use. When not given or empty, the Temperature
      %      is searched in the object. The temperature is in [K]. 1 meV=11.605 K.
      % output:
      %   sqw: S(q,w) 2D data set (iData_Sqw2D)
      %
      % conventions:
      % w = omega = Ei-Ef = energy lost by the neutron [meV]
      %    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
      %    omega < 0, neutron gains energy, anti-Stokes
      %
      % references: M. Mattes and J. Keinert, IAEA INDC(NDS)-0470, 2005.
      %             R. E. MacFarlane, LA-12639-MS (ENDF-356), 1994.
      %
      % Example: sab = iData_Sab('SQW_coh_lGe.nc');
      %          sqw = Sqw(sab,72.6,1235);
      %          sab2= Sab(sqw);
      %          subplot(log([Sqw Sab Sqw2))
      %
      % (c) E.Farhi, ILL. License: EUPL.
      s = Sab2Sqw(self);  % private
    end
    
    function f = iData(self)
      % iData_Sab: iData(sab): convert a iData_Sab back to iData
      %   The specific iData_Sab methods are then not accessible anymore.
      %
      % syntax:
      %   d = iData(sab)
      %
      % input:  sab: S(alpha,beta) [iData_Sab]
      % output: iData
      %
      % sab=iData_Sab('SQW_coh_lGe.nc', 72.6,1235);  iData(sab)
      %
      % See also: iData, iData_Sqw2D
      f = [];
      for index=1:numel(self)
        this = self(index);
        f1   = iData;
        this = struct(this);
        w = warning('query','iData:subsasgn');
        warning('off','iData:subsasgn');
        for p = fieldnames(this)'
          f1.(p{1}) = this.(p{1}); % generates a warning when setting Alias
        end
        warning(w);
        f = [ f f1 ];
      end
    end
    
    function DOS = dos(self, varargin)
      % iData_Sab: dos(sab): compute the generalised density of states (gDOS)
      %  The gDOS is an approximation of the vibrational spectra (DOS).
      %  This routine should better be applied on an incoherent dynamic S(a,b) data set.
      %  The incoherent approximation states that the gDOS from an incoherent S(a,b) is 
      %  roughly equal to that obtained from a coherent S(a,b).
      %
      %       gDOS(q,w) = S(q,w) w^2/q^2                   Bellissent
      %       gDOS(q,w) = S(q,w) w  /q^2/[1 + n(hw)]       Carpenter/Price
      %  and:
      %       gDOS(w)   = lim(q->0) [ gDOS(q,w) ]
      %
      %       gDOS(q,w) = w*S(q,w)/(1+n(w))/(Q^4max-Q^4min) Bredov/Oskotskii
      %       gDOS(w)   = trapz(g, 2)
      %
      %  The Bredov/Oskotskii methodology provides the best gDOS estimate, using the
      %  whole data set.
      %
      %  The applicability to a coherent dynamic structure factor S(q,w) should be
      %  taken with great care, as this formalism then does not fully hold.
      %
      %  The method to use in the gDOS computation can be given as 2nd argument
      %       gDOS = dos(Sqw, 'Bredov')         more accurate as it uses 100% of data
      %       gDOS = dos(Sqw, 'Carpenter')      Temperature must be a property
      %       gDOS = dos(Sqw, 'Bellissent')     simple yet efficient
      %
      %  The gDOS is stored in the 'gDOS' property, and is retrieved without recomputation
      %  when available. To force a recomputation of the gDOS, use:
      %       dos(Sab, 'method', 0) or dos(Sab, 'method', 'force')
      %
      %   When no output is used, a plot is shown.
      %
      % References: Price J. et al, Non Cryst Sol 92 (1987) 153
      %         Bellissent-Funel et al, J. Mol. Struct. 250 (1991) 213
      %         Carpenter and Pelizarri, Phys. Rev. B 12, 2391 (1975)
      %         Suck et al, Journal of Alloys and Compounds 342 (2002) 314
      %         Bredov et al., Sov. Phys. Solid State 9, 214 (1967)
      %         V.S. Oskotskii, Sov. Phys. Solid State 9 (1967), 420.
      %
      % syntax:
      %   g = dos(sab)
      %   g = dos(sab, method, n)
      %
      % input:  sab:    S(alpha,beta) [iData_Sab]
      %         method: 'Carpenter','Bellissent' or 'Bredov' (default)
      %         n:      number of low-angle values to integrate (integer). Default is 10.
      %                 can also be 0 or 'force' to re-compute the gDOS.
      % output: g:      gDOS [iData_vDOS]
      %
      % Example: sab=iData_Sab('SQW_coh_lGe.nc',72.6,1235); dos(sab);
      %
      % See also: iData_Sqw2D/dos, iData_vDOS
      DOS = dos(iData_Sqw2D(self), varargin{:});
      setalias(self, 'gDOS', DOS, DOS.Title);
      if length(inputname(1))
        assignin('caller',inputname(1),self);
      end
      
      if nargout == 0 && ~isempty(DOS)
        fig=figure; 
        % plot total DOS
        h=plot(DOS); set(h,'LineWidth',2);
        set(fig, 'NextPlot','new');
      end
    end
    
    function t = thermochemistry(s, T, options)
      % iData_Sab: thermochemistry(sab, T): compute thermodynamic quantities for 2D S(alpha,beta) data sets.
      %   When no output is used, a plot is shown.
      %
      % The function returns an array of iData objects:
      %   density of states
      %   entropy                           S [eV/K/cell]
      %   internal_energy                   U [eV/cell]
      %   helmholtz_energy                  F [eV/cell]
      %   heat_capacity at constant volume Cv [eV/K/cell]
      %
      % When no output is used, a plot is shown.
      %
      % syntax:
      %   t = thermochemistry(sab)
      %   t = thermochemistry(sab, Tmin:Tmax)
      %
      % input:  sab:  S(alpha,beta) [iData_Sab]
      %         T:    Temperature [scalar or vector]. Default 1-500.
      % output: DOS,S,U,F,Cv [iData_array]
      %
      % Example: sab=iData_Sab('SQW_coh_lGe.nc',72.6,1235); thermochemistry(sab);
      %
      % See also: iData_Sab/dos, iData_vDOS/thermochemistry
      if nargin < 2, T=[]; end
      if nargin < 3, options=[]; end
      if isempty(T)
        T = 1:500;  % default
      end
      
      if nargout == 0 || ~isempty(strfind(options, 'plot'))
        options = [ options ' newplot' ];
      end
      s = iData_Sqw2D(s);
      t = thermochemistry(s, T, options);
    end
    
    function sigma = moments(data, varargin)
      % iData_Sab: moments(sab, M, T, classical): compute Sqw moments/sum rules (harmonic frequencies)
      %
      % Compute the structure factor (moment 0), recoil energy (moment 1) and the
      %   collective, harmonic and mean energy transfer dispersions.
      %
      % The result is given as an iData array with data sets:
      %   S(q) = \int S(q,w) dw = <S(q,w)>                 structure factor [moment 0]
      %   Er   = \int w*S(q,w) dw = <wS(q,w)> = h2q2/2M       recoil energy [moment 1]
      %   Wc   = sqrt(2kT*Er/S(q))                    collective/isothermal dispersion
      %   Wl                                          harmonic/longitudinal excitation
      %   Wq   = 2q*sqrt(kT/S(q)/M)                               mean energy transfer
      %   M2   = <w2S(q,w)>                                                 [moment 2]
      %   M3   = <w3S(q,w)>                                                 [moment 3]
      %   M4   = <w4S(q,w)>                                                 [moment 4]
      %
      % When no output is used, a plot is shown.
      %
      % Reference: 
      %   Helmut Schober, Journal of Neutron Research 17 (2014) pp. 109
      %   Lovesey, Theory of Neutron Scattering from Condensed Matter, Vol 1, p180 eq. 5.38 (w0)
      %   J-P.Hansen and I.R.McDonald, Theory of simple liquids Academic Press New York 2006.
      %
      % syntax:
      %   m = moments(sab)
      %   m = moments(sab, M, T, classical)
      %
      % input: sab: S(alpha,beta) [iData_Sab]
      %        M:   Molar weight [g/mol]. Searched in the input data set when not given.
      %        T:   Temperature [K]. Searched in the input data set when not given.
      %        classical: 0 for non symmetric S(a,b) [with Bose, from exp.], 1 for symmetric (from MD)
      % output:
      %   moments=[ sq M1 wc wl wq M2 M3 M4 ] as an iData array
      %
      % Example: sab=iData_Sab('SQW_coh_lGe.nc',72.6,1235); moments(sab);
      %
      % See also: iData_Sqw2D/moments
      sigma = moments(iData_Sqw2D(data), varargin{:});
      if nargout == 0
        fig=figure; h  =plot(sigma(1:6),'legend'); set(fig, 'NextPlot','new');
      end
    end
    
    function s = symmetrize(s0)
      % iData_Sab: symmetrize(s): extend the S(alpha,beta) in both beta/energy sides
      %  The resulting S(alpha,beta) is the combination of S(alpha,beta) and 
      %  S(alpha,-beta), which is thus symmetric in beta(energy):
      %     S(alpha,beta) = S(alpha,-beta)
      %
      %  The incoming data set should NOT contain the Bose factor, that is it
      %    should be 'classical'.
      %  To obtain a 'classical' S(alpha,beta) from an experiment, use first:
      %    deBosify(s)
      %
      % The positive energy values in the S(alpha,-beta) map correspond to Stokes 
      %   processes i.e. material gains energy, and neutrons loose energy when 
      %   down-scattered. Neutron up-scattering corresponds to anti-Stokes processes.
      %
      % When no output is used, a plot is shown.
      %
      % syntax:
      %   sab_sym = symmetrize(sab)
      %
      % input:
      %   s:  Sab data set (classical, often labelled as S*)
      % output:
      %   s:  S(alpha,bela) symmetrised in energy, classical.
      %
      % Example: sab=iData_Sab('SQW_coh_lGe.nc',72.6,1235); symmetrize(sab);
      %
      % See also: iData_Sab/Bosify, iData_Sab/deBosify, 
      %           iData_Sab/dynamic_range, iData_Sab/scattering_cross_section
      s = Sab(symmetrize(Sqw(s0)));
      s.Title = [ 'symmetrize(' s0.Title  ')' ];
      title(s, [  'symmetrize(' title(s0) ')' ]);
      s.Label = [ 'symmetrize(' s0.Label  ')' ];
      if nargout == 0
        fig=figure; h  =plot(log10(s)); set(fig, 'NextPlot','new');
      end
    end
    
    function s = Bosify(s0, varargin)
      % iData_Sab: Bosify(s, T): apply the 'Bose' factor (detailed balance) to a classical data set.
      %   The initial data set should obey S*=S(a,b) = S(a,-b), i.e. be 'classical'.
      %   The resulting data set is 'quantum/experimental' and satisfies the detailed 
      %   balance. It contains the temperature effect (population).
      %
      % conventions:
      % beta = (Ef-Ei)/kT = -hw/kT  energy gained by the neutron, unitless
      %    beta < 0, neutron looses energy, down-scattering (Stokes)
      %    beta > 0, neutron gains energy,  up-scattering (anti-Stokes)
      %
      %    S(q,-w) = exp(-hw/kT) S(q,w)
      %    S(q,w)  = exp( hw/kT) S(q,-w)
      %    S(q,w)  = Q(w) S*(q,w) with S*=classical limit and Q(w) defined below.
      % For beta > 0, S(q,-beta) > S(q,beta)
      %               
      % The semi-classical correction, Q, aka 'quantum' correction factor, 
      % can be selected from the optional   'type' argument:
      %    Q = exp(hw_kT/2)                 'Schofield' or 'Boltzmann'
      %    Q = hw_kT./(1-exp(-hw_kT))       'harmonic'  or 'Bader'
      %    Q = 2./(1+exp(-hw_kT))           'standard'  or 'Frommhold' (default)
      %
      % The 'Boltzmann' correction leads to a divergence of the S(q,w) for e.g. w above 
      % few 100 meV. The 'harmonic' correction provides a reasonable correction but does
      % not fully avoid the divergence at large energies.
      %
      %  Bose factor: n(w) = 1./(exp(w*11.605/T) -1) ~ exp(-w*11.605/T)
      %               w in [meV], T in [K]
      %
      % References:
      %  B. Hehr, http://www.lib.ncsu.edu/resolver/1840.16/7422 PhD manuscript (2010).
      %  S. A. Egorov, K. F. Everitt and J. L. Skinner. J. Phys. Chem., 103, 9494 (1999).
      %  P. Schofield. Phys. Rev. Lett., 4, 239 (1960).
      %  J. S. Bader and B. J. Berne. J. Chem. Phys., 100, 8359 (1994).
      %  T. D. Hone and G. A. Voth. J. Chem. Phys., 121, 6412 (2004).
      %  L. Frommhold. Collision-induced absorption in gases, 1 st ed., Cambridge
      %    Monographs on Atomic, Molecular, and Chemical Physics, Vol. 2,
      %    Cambridge Univ. Press: London (1993).
      %
      % syntax:
      %   sab_T = Bosify(sab)
      %
      % input: s:     Sab data set (classical, symmetric in energy, no T Bose factor)
      %        T:     Temperature [K]. Searched in the input data set when not given.
      %        type:  'Schofield' or 'harmonic' or 'standard' (default)
      %
      % output:
      %   sb: quantum Sab data set (experimental/non classical, iData_Sab).
      %
      % Example: sab=iData_Sab('SQW_coh_lGe.nc',72.6,1235); sb=Bosify(symmetrize(sab));
      %
      % See also: iData_Sab/deBosify, iData_Sab/symmetrize, iData_Sqw2D
      s = iData_Sab(Bosify(iData_Sqw2D(s0), varargin{:}));
      s.Title = [ 'Bosify(' s0.Title  ')' ];
      title(s, [  'Bosify(' title(s0) ')' ]);
      s.Label = [ 'Bosify(' s0.Label  ')' ];
      setalias(s,'classical',   0, '[0=from measurement, with Bose factor included, 1=from MD, symmetric]');
      if nargout == 0
        fig=figure; h  =plot(log10(s)); set(fig, 'NextPlot','new');
      end
    end
    
    function s = deBosify(s0, varargin)
      % iData_Sab: deBosify(s, T): remove the 'Bose' factor (detailed balance) from an experimental/ quantum data set.
      %   The initial data set is 'quantum/experimental' and satisfies the detailed balance.
      %   The resulting data set obeys S*=S(a,b) = S(a,-b), i.e. is 'classical'. 
      %   It suppresses the temperature effect (population).
      %
      % conventions:
      % beta = (Ef-Ei)/kT = -hw/kT  energy gained by the neutron, unitless
      %    beta < 0, neutron looses energy, down-scattering (Stokes)
      %    beta > 0, neutron gains energy,  up-scattering (anti-Stokes)
      %
      %    S(q,-w) = exp(-hw/kT) S(q,w)
      %    S(q,w)  = exp( hw/kT) S(q,-w)
      %    S(q,w)  = Q(w) S*(q,w) with S*=classical limit and Q(w) defined below.
      % For beta > 0, S(q,-beta) > S(q,beta)
      %               
      % The semi-classical correction, Q, aka 'quantum' correction factor, 
      % can be selected from the optional   'type' argument:
      %    Q = exp(hw_kT/2)                 'Schofield' or 'Boltzmann'
      %    Q = hw_kT./(1-exp(-hw_kT))       'harmonic'  or 'Bader'
      %    Q = 2./(1+exp(-hw_kT))           'standard'  or 'Frommhold' (default)
      %
      % The 'Boltzmann' correction leads to a divergence of the S(q,w) for e.g. w above 
      % few 100 meV. The 'harmonic' correction provides a reasonable correction but does
      % not fully avoid the divergence at large energies.
      %
      %  Bose factor: n(w) = 1./(exp(w*11.605/T) -1) ~ exp(-w*11.605/T)
      %               w in [meV], T in [K]
      %
      % References:
      %  B. Hehr, http://www.lib.ncsu.edu/resolver/1840.16/7422 PhD manuscript (2010).
      %  S. A. Egorov, K. F. Everitt and J. L. Skinner. J. Phys. Chem., 103, 9494 (1999).
      %  P. Schofield. Phys. Rev. Lett., 4, 239 (1960).
      %  J. S. Bader and B. J. Berne. J. Chem. Phys., 100, 8359 (1994).
      %  T. D. Hone and G. A. Voth. J. Chem. Phys., 121, 6412 (2004).
      %  L. Frommhold. Collision-induced absorption in gases, 1 st ed., Cambridge
      %    Monographs on Atomic, Molecular, and Chemical Physics, Vol. 2,
      %    Cambridge Univ. Press: London (1993).
      %
      % syntax:
      %   sab = deBosify(sab_T)
      %
      % input: sab_T: Sab data set (quantum/experimental, with T Bose factor)
      %        T:     Temperature [K]. Searched in the input data set when not given.
      %        type:  'Schofield' or 'harmonic' or 'standard' (default)
      %
      % output:
      %   sb: classical/symmetric Sab data set (iData_Sab).
      %
      % Example: sab=iData_Sab('SQW_coh_lGe.nc',72.6,1235); sb=Bosify(symmetrize(sab));
      %          sab0=deBosify(sb);
      %
      % See also: iData_Sab/Bosify, iData_Sab/symmetrize, iData_Sqw2D
      s = iData_Sab(deBosify(iData_Sqw2D(s0), varargin{:}));
      s.Title = [ 'deBosify(' s0.Title  ')' ];
      title(s, [  'deBosify(' title(s0) ')' ]);
      s.Label = [ 'deBosify(' s0.Label  ')' ];
      setalias(s,'classical',   1, '[0=from measurement, with Bose factor included, 1=from MD, symmetric]');
      if nargout == 0
        fig=figure; h  =plot(log10(s)); set(fig, 'NextPlot','new');
      end
    end
    
    function [s, sphiw] = dynamic_range(s, varargin)
      % iData_Sab: dynamic_range: crop the S(alpha,beta) to the available dynamic range
      %   for given incident neutron energy.
      %
      % The dynamic range is defined from the momentum and energy conservation laws:
      %  Ef         = Ei - w                                is positive
      %  cos(theta) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf)  is within [-1:1]
      %
      % The incident neutron energy can be computed using:
      %  Ei = 2.0721*Ki^2 = 81.8042/lambda^2 with Ki in [Angs-1] and lambda in [Angs]
      %
      % conventions:
      % beta = (Ef-Ei)/kT = -hw/kT  energy gained by the neutron, unitless
      %    beta < 0, neutron looses energy, down-scattering (Stokes)
      %    beta > 0, neutron gains energy,  up-scattering (anti-Stokes)
      %
      % The scattering angle phi can be restricted to match a detection area
      % with the syntax:
      %   sab_Ei=dynamic_range(sab, Ei, [angles])
      %
      % The syntax:
      %   [sab_Ei, sphiw] = dynamic_range(sab,...)
      % also returns the S(phi,w) data set, which shows the detected signal vs scattering 
      % angle and energy transfer. This is the double differential cross section 
      % as a function of the scattering angle.
      %
      % syntax:
      %   sab_Ei        =dynamic_range(sab, Ei)
      %   [sab_Ei,sphiw]=dynamic_range(sab, Ei)
      %   [sab_Ei,sphiw]=dynamic_range(sab, Ei, [min_angle max_angle])
      %
      % input:
      %   sab:    S(alpha,beta) data set
      %   Ei:     incoming neutron energy [meV]
      %   angles: detection range in [deg] as a vector. Min and Max values are used.
      % output:
      %   sab_Ei: S(alpha,beta) cropped to dynamic range for incident energy Ei.
      %   sphiw:  S(phi,w) angular dynamic structure factor for incident energy Ei.
      %
      % Example: sab=iData_Sab('SQW_coh_lGe.nc',72.6,1235); dynamic_range(symmetrize(sab), 14.8, [-20 135])
      %
      % See also: iData_Sab/Bosify, iData_Sab/deBosify, iData_Sab/symmetrize, iData_Sab/scattering_cross_section

      [s, sphiw] = dynamic_range(iData_Sqw2D(s), varargin{:});
      s = iData_Sab(s);
      if nargout == 0
        fig=figure; h  =subplot(log10([ s sphiw ]), 'view2 tight'); set(fig, 'NextPlot','new');
      end
    end
    
    function sigma = scattering_cross_section(s, varargin)
      % iData_Sab: scattering_cross_section: compute the total neutron scattering cross section for
      %   incoming neutron energy. The S(alpha,beta) should be the non-classical
      %   dynamic structure factor. 
      %
      %   Such data sets are obtained from e.g. xray and neutron scattering 
      %   experiments on isotropic density materials (liquids, powders, amorphous
      %   systems). 
      %
      %   Data sets from analytical models and molecular dynamics simulations must 
      %   be symmetrised in energy, and the detailed balance must be applied to 
      %   take into account the material temperature on the inelastic part.
      %
      %   The incident neutron energy is given in [meV], and may be computed:
      %     Ei = 2.0721*Ki^2 = 81.8042/lambda^2 with Ki in [Angs-1] and lambda in [Angs]
      %     
      %   The S(alpha,beta) is first restricted to the achievable dynamic range:
      %     Ef         = Ei - w                                is positive
      %     cos(theta) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf)  is within [-1:1]
      %   and then integrated as XS = 1/2Ki^2 \int q S(q,w) dq dw
      %
      %   The computed value must then be multiplied by the tabulated bound cross-section
      %   from e.g. Sears, Neut. News 3 (1992) 26.
      %
      %   When the weight M of the scattering unit is given, it is used to multiply the
      %    cross section by the estimated Debye-Waller-like factor so that it equals
      %    A/(A+1)]^2 at Ei=1eV to gradually go from the bound (thermal) to the free 
      %    cross section at 1 eV. The threshold of 1eV is used by e.g. OpenMC.
      %     W  = 2e-3*(log(M)-log(M+1));
      %     DW = exp(W*Ei) = exp(2e-3*(log(M)-log(M+1))*Ei)
      %   Above the epithermal energy threshold, the DW factor is kept fixed. 
      %   For a poly-atomic scatterer, the effective mass is computed by weighting
      %   with the bound cross sections for elements A, e.g.
      %     r = sqrt(sum((A/(A+1))^2 * sigma_bound(A)) / sum(sigma_bound(A)));
      %     M = r/(1-r)
      %   For instance, for H2O (twice A=1 sigma=80.2 and A=18 sigma=4.2):
      %     r = sqrt((2*(1/(1+1))^2*80.2+(18/(18+1))^2*4.2)/(2*80.2+4.2)) = 0.52
      %     M = r/(1-r) = 1.06 i.e. scattering is mostly the hydrogen one.
      %   WARNING: this factor is NOT the Debye-Waller factor exp(-<u2>Q2) !
      %
      % A classical S(a,b) obeys S(a,b) = S(a,-b).
      % The non classical S(a,b), needed by this function, can be obtained from a 
      % classical S(a,b) (which is symmetric in beta) with e.g.:
      %   extend to +/- energy range
      %     s = symmetrize(s); 
      %   apply detailed balance (Bose factor). Omit T if you do not know it.
      %     s = Bosify(s, T);
      %
      % The positive energy values in the S(q,w) map correspond to Stokes processes, 
      % i.e. material gains energy, and neutrons loose energy when down-scattered.
      %
      % syntax:
      %   sigma = scattering_cross_section(sab, Ei)
      %   sigma = scattering_cross_section(sab, Ei, M)
      %
      % input:
      %   sab: S(alpha,beta) data set (non classical, with T Bose factor e.g from experiment)
      %        e.g. 2D data set with w as 1st axis (rows, meV), q as 2nd axis (Angs-1).
      %   Ei: incoming neutron energy [meV]
      %   M: molar weight of the atom/molecule in [g/mol].
      %     when given empty, it is searched 'weight' or 'mass' is the object.
      %     Default is set to 0, i.e. the Debye-Waller factor is not taken into account.
      % output:
      %   sigma: cross section per scattering unit (scalar or iData)
      %          to be multiplied afterwards by the bound cross section [barn]
      %
      % Example: sab=iData_Sab('SQW_coh_lGe.nc',72.6,1235); 
      %          sigma = scattering_cross_section(Bosify(symmetrize(sab)), 14.6);
      %
      % See also: iData_Sab/Bosify, iData_Sab/deBosify, iData_Sab/symmetrize, iData_Sab/dynamic_range
      sigma = scattering_cross_section(iData_Sqw2D(s), varargin{:});
      if nargout == 0
        fig=figure; h  =plot(sigma); set(fig, 'NextPlot','new');
      end
    end
    
    function sq = structure_factor(s)
      % iData_Sab: structure_factor: compute the structure factor
      %  The structure factor is the integral of the dynamic structure factor along
      %  the energy axis. It is representative of the material structure.
      %  Its Fourier transform is the pair distribution function g(r).
      %
      %  This function is basically a call to trapz(s,1)
      %
      % References: Fischer, Barnes and Salmon, Rev Prog Phys 69 (2006) 233
      %
      % syntax:
      %   sq = structure_factor(s)
      %
      % input:
      %   s:  S(alpha,beta) data set
      % output:
      %   sq: S(q) data set
      %
      % Example: sq = structure_factor(iData_Sab('SQW_coh_lGe.nc',72.6,1235));
      %
      % See also: trapz, iData/trapz, iData_Sqw2D/structure_factor
      sq = structure_factor(iData_Sqw2D(s));
      if nargout == 0
        fig=figure; h  =plot(sq); set(fig, 'NextPlot','new');
      end
    end
    
    function [inc, multi] = incoherent(s, varargin)
      % iData_Sab: multi_phonons_incoherent: compute the multi-phonon contributions in S(q,w) from an initial density of states in the incoherent gaussian approximation
      %
      % This implementation is in principle exact for an isotropic monoatomic material,
      % e.g. a liquid or powder.
      % This methodology is equivalent to the LEAPR module of NJOY ("phonon expansion")
      % to compute S(alpha,beta) from a vibrational density of states.
      %
      % Reference:
      %   H. Schober, Journal of Neutron Research 17 (2014) 109–357
      %     DOI 10.3233/JNR-140016 (see esp. pages 328-331)
      %   V.S. Oskotskii, Sov. Phys. Solid State 9 (1967), 420.
      %   A. Sjolander, Arkiv for Fysik 14 (1958), 315.
      %
      % syntax:
      %   [Sinc, multi] = multi_phonons_incoherent(gw)
      %   [Sinc, multi] = multi_phonons_incoherent(gw, q, T, sigma, m, n)
      %
      % input:
      %   gw: the vibrational density of states per [meV] [iData]
      %   q:  the momentum axis [Angs-1, vector]
      %   T:  temperature [K]
      %   sigma: neutron cross section [barns]
      %   m:  mass [g/mol]
      %   n:  number of iterations in the series expansion, e.g. 5
      %
      % output:
      %   Sinc:   neutron weighted incoherent S(q,w) single-phonon term [iData_Sqw2D]
      %   multi:  so called phonon expansion aka multi-phonon terms
      %
      % Example:
      %   sab=iData_Sab('SQW_coh_lGe.nc',72.6,1235); incoherent(sab)
      %
      % See also: iData_Sab/dos, iData_Sab/multi_phonons
      [inc,multi] = incoherent(Sqw(s), varargin{:});
      inc = Sab(inc);
      multi=Sab(multi);
      if nargout == 0
        fig=figure; 
        h  =subplot(log10([inc multi]), 'view2'); 
        set(fig, 'NextPlot','new');
      end
    end
    
    function [G1, G2] = multi_phonons(s, varargin)
      % iData_Sqw2D: multi_phonons: compute the integrated multi-phonon DOS from an initial density of states
      %
      % The 'generalized' neutron weighted density of states (gDOS) is computed from an 
      % initial vibrational density of states (vDOS).
      % Missing arguments (or given as [] empty), are searched in the initial density 
      % of states.
      %
      % This implementation is in principle exact for an isotropic monoatomic material,
      % e.g. a liquid or powder.
      % This methodology is a complete rewrite of the MUPHOCOR code.
      %
      % Reference:
      %   H. Schober, Journal of Neutron Research 17 (2014) 109–357
      %     DOI 10.3233/JNR-140016 (see esp. pages 328-331)
      %   V.S. Oskotskii, Sov. Phys. Solid State 9 (1967), 420.
      %   A. Sjolander, Arkiv for Fysik 14 (1958), 315.
      %   W. Reichardt, MUPHOCOR Karlsruhe Report 13.03.01p06L (1984)
      %
      % syntax:
      %   [G,multi] = multi_phonons(sab)
      %   [G,multi] = multi_phonons(sab, Ki, T, sigma, m, n)
      %
      % input:
      %   sab:  S(alpha,beta) data set
      %   Ki:   incident wavevector [Angs-1]. Ki=2*pi/lambda=0.695*sqrt(Ei)
      %   T:    temperature [K]. Searched in the input data set when not given.
      %   sigma: neutron cross section [barns]. Searched in the input data set when not given.
      %   m:    mass [g/mol]. Searched in the input data set when not given.
      %   phi:  detector angles, min and max [deg]
      %   n:    number of iterations in the series expansion, e.g. 5
      %
      % output:
      %   G:      total gDOS [p=1...]
      %   multi:  multi-phonon gDOS contribution [p=2...]
      %
      % Example: sab=iData_Sab('SQW_coh_lGe.nc',72.6,1235); multi_phonons(sab)
      %
      % See also: iData_Sab/dos, iData_Sab/incoherent
      
      g  = dos(s);
      G1 = multi_phonons_dos(g, varargin{:});
      G2 = plus(G1(2:end)); title(G2, 'gDOS [multi]');        G2.Label=title(G2); G2.Title=title(G2);
      G1 = plus(G1);        title(G1, 'gDOS [single+multi]'); G1.Label=title(G1); G1.Title=title(G1);
      if nargout == 0
        fig=figure; 
        h  =plot([ G1 G2 ]); 
        set(fig, 'NextPlot','new');
      end
    end
    
  end % methods
  
end
