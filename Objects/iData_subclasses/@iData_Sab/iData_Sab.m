classdef iData_Sab < iData
  % iData_Sab: create a 2D S(alpha,beta) data set (iData flavour)
  %
  % The iData_Sab class is a 2D data set holding a S(alpha,beta) dynamic 
  %   structure factor aka scattering function/law.
  %   The data set axes are beta as 1st axis (rows), alpha as 2nd axis (Angs-1).
  %
  %  The S(alpha,beta) is a representation of the dynamic structure factor 
  %  using unitless momentum and energy variables defined as:
  %     alpha= h2q2/2MkT  = (Ei+Ef-2*mu*sqrt(Ei*Ef))/AkT
  %     beta = -hw/kT     = (Ef-Ei)/kT                  energy gained by neutron
  %     A    = M/m
  %     mu   = cos(theta) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf)
  %  This representation is common in nuclear data, neutron sections (e.g. ENDF MF7)
  %
  % Example: sab=iData_Sab('SQW_coh_lGe.nc', 72.6,300)
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
  % inc = incoherent(sab)
  %   Compute the incoherent neutron scattering law estimate in the incoherent 
  %     gaussian approximation, using its density of states.
  %
  % [G,M] = multi_phonons(sab)
  %   Compute the integrated multi-phonon DOS terms from an initial density of states
  %
  % input:
  %   can be an iData or filename to generate a 2D Sab object.
  %
  % output: an iData_Sab object
  %
  % See also: iData
  % (c) E.Farhi, ILL. License: EUPL.
  
  properties
  end
  
  methods
    % main instantiation method
    function obj = iData_Sab(s, varargin)
      % iData_Sab: create the iData_Sab subclass
      %
      % input:
      %   s: must be a 2D iData object holding a S(q,w), S(phi,t), S(alpha,beta) or S(phi,w)
      %      or file name with such data.
      %   M: molar mass [g/mol]. Searched in the input data set when not given.
      %   T: temperature [K]. Searched in the input data set when not given.
      % output:
      %   S(alpha,beta) data set
      %
      % Example: sab=iData_Sab('SQW_coh_lGe.nc', 72.6,300);
      %
      % See also: iData_Sqw2D
      obj = obj@iData;
      obj.class = mfilename;
      
      if ~nargin, return; end  % empty object
      
      % convert/test
      if isa(s, 'iData_Sqw2D') s = Sab(s); end
      if ~isa(s, 'iData'), m = Sab(iData_Sqw2D(s), varargin{:}); else m=s; end
      if ~isa(m, 'iData') || isempty(m) || ndims(m) ~= 2
        error([ mfilename ': the given input ' class(s) ' does not seem to be convertible to iData_Sab.' ])
      end
      
      % transfer properties
      % this is a safe way to instantiate a subclass
      warning off MATLAB:structOnObject
      m = struct(m);
      for p = fieldnames(m)'
        obj.(p{1}) = m.(p{1});
      end
      obj.class = mfilename;
 
    end % iData_Sab constructor
    
    % parameters (search for parameters in iData)
    function parameters = parseparams(s)
      % iData_Sab: parseparams(sab): search for physical quantities in object.
      %
      % input:  sab: S(alpha,beta) [iData_Sab]
      % output: physical parameters [struct]
      %
      % See also: iData_Sqw2D/parseparams
      [s,parameters,fields] = Sqw_parameters(s);
      if nargout == 0 & length(inputname(1))
        assignin('caller',inputname(1),s);
      end
    end
    
    function f = iData(self)
      % iData_Sab: iData(sab): convert a iData_Sab back to iData
      %
      % input:  sab: S(alpha,beta) [iData_Sab]
      % output: iData
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
      %   The gDOS is an approximation of the vibrational spectra (DOS).
      %   When no output is used, a plot is shown.
      %
      % input:  sab:    S(alpha,beta) [iData_Sab]
      %         method: 'Carpenter','Bellissent' (default) or 'Bredov'
      % output: gDOS [iData_vDOS]
      %
      % Example: dos(iData_Sab('SQW_coh_lGe.nc',72.6,300));
      %
      % See also: iData_Sqw2D/dos, iData_vDOS
      DOS = dos(iData_Sqw2D(self));
      
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
      % input:  sab:  S(alpha,beta) [iData_Sab]
      %         T:    Temperature [scalar or vector]. Default 1-500.
      % output: DOS,S,U,F,Cv [iData_array]
      %
      % Example: thermochemistry(iData_Sab('SQW_coh_lGe.nc',72.6,300));
      %
      % See also: iData_Sab/dos, iData_vDOS/thermochemistry
      if nargin < 2, T=[]; end
      if nargin < 3, options=[]; end
      if isempty(T)
        T = 1:500;  % default
      end
      
      if nargout == 0 || ~isempty(strfind(options, 'plot'))
        options = [ options ' plot' ];
      end
      s = iData_Sqw2D(s);
      t = thermochemistry(s, T, options);
    end
    
    function sigma = moments(data, varargin)
      % iData_Sab: moments(sab, M, T, classical): compute Sqw moments/sum rules (harmonic frequencies)
      %   When no output is used, a plot is shown.
      %
      % input: sab: S(alpha,beta) [iData_Sab]
      %        M:   Molar weight [g/mol]. Searched in the input data set when not given.
      %        T:   Temperature [K]. Searched in the input data set when not given.
      %        classical: 0 for non symmetric S(a,b) [with Bose, from exp.], 1 for symmetric (from MD)
      % output:
      %   moments=[ sq M1 wc wl wq M2 M3 M4 ] as an iData array
      %
      % Example: moments(iData_Sab('SQW_coh_lGe.nc',72.6,300));
      %
      % See also: iData_Sqw2D/moments
      sigma = moments(iData_Sqw2D(data), varargin{:});
      if nargout == 0
        fig=figure; h  =plot(sigma(1:6),'legend'); set(fig, 'NextPlot','new');
      end
    end
    
    function s = symmetrize(s)
      % iData_Sab: symmetrize(s): extend the S(alpha,beta) in both beta/energy sides
      %  The resulting S(alpha,beta) is the combination of S(alpha,beta) and 
      %  S(alpha,-beta), which is thus symmetric in beta(energy):
      %     S(alpha,beta) = S(alpha,-beta)
      %
      %  The incoming data set should NOT contain the Bose factor, that is it
      %    should be 'classical'.
      %  To obtain a 'classical' S(alpha,beta) from an experiment, use first:
      %    deBosify(s, T)
      %
      % The positive energy values in the S(alpha,-beta) map correspond to Stokes , 
      % processes i.e. material gains energy, and neutrons loose energy when 
      % down-scattered. Neutron up-scattering corresponds to anti-Stokes processes.
      %
      %   When no output is used, a plot is shown.
      %
      % input:
      %   s:  Sab data set (classical, often labelled as S*)
      % output:
      %   s:  S(alpha,bela) symmetrised in energy, classical.
      %
      % Example: sab=iData_Sab('SQW_coh_lGe.nc',72.6,300); symmetrize(sab);
      %
      % See also: iData_Sab/Bosify, iData_Sab/deBosify, 
      %           iData_Sab/dynamic_range, iData_Sab/scattering_cross_section
      s = Sab(symmetrize(Sqw(s)));
      if nargout == 0
        fig=figure; h  =plot(s); set(fig, 'NextPlot','new');
      end
    end
    
    function s = Bosify(s0, varargin)
      % iData_Sab: Bosify(s, T): apply the 'Bose' factor (detailed balance) to a classical data set.
      %   The initial data set should obey S*=S(a,b) = S(a,-b), i.e. be 'classical'.
      %
      % input: s:     Sab data set (classical, symmetric in energy, no T Bose factor)
      %        T:     Temperature [K]. Searched in the input data set when not given.
      %        type:  'Schofield' or 'harmonic' or 'standard' (default)
      %
      % conventions:
      % beta = (Ef-Ei)/kT = -hw/kT  energy gained by the neutron, unitless
      %    beta < 0, neutron looses energy, down-scattering (Stokes)
      %    beta > 0, neutron gains energy,  up-scattering (anti-Stokes)
      % Egelstaff, Eq (9.25) p189,  with Q(w) is defined below:
      %    S(q,-w) = exp(-hw/kT) S(q,w)
      %    S(q,w)  = exp( hw/kT) S(q,-w)
      %    S(q,w)  = Q(w) S*(q,w) with S*=classical limit
      % for omega > 0, S(q,w) > S(q,-w)
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
      % Example: sab=iData_Sab('SQW_coh_lGe.nc',72.6,300); sb=Bosify(symmetrize(sab));
      %
      % See also: iData_Sab/deBosify, iData_Sab/symmetrize, iData_Sqw2D
      s = iData_Sab(Bosify(iData_Sqw2D(s0), varargin{:}));
      setalias(s,'classical',   0, '[0=from measurement, with Bose factor included, 1=from MD, symmetric]');
    end
    
    function s = deBosify(s0, varargin)
      % iData_Sab: deBosify(s, T): remove Bose factor (detailed balance) from an 'experimental/quantum' data set.
      %   The resulting data set obeys S*=S(a,b) = S(a,-b), i.e. is 'classical'.
      %
      % input: s:     Sab data set (quantum/experimental, with T Bose factor)
      %        T:     Temperature [K]. Searched in the input data set when not given.
      %        type:  'Schofield' or 'harmonic' or 'standard' (default)
      %
      % conventions:
      % beta = (Ef-Ei)/kT = -hw/kT  energy gained by the neutron, unitless
      %    beta < 0, neutron looses energy, down-scattering (Stokes)
      %    beta > 0, neutron gains energy,  up-scattering (anti-Stokes)
      % Egelstaff, Eq (9.25) p189,  with Q(w) is defined below:
      %    S(q,-w) = exp(-hw/kT) S(q,w)
      %    S(q,w)  = exp( hw/kT) S(q,-w)
      %    S(q,w)  = Q(w) S*(q,w) with S*=classical limit
      % for omega > 0, S(q,w) > S(q,-w)
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
      % Example: sab=iData_Sab('SQW_coh_lGe.nc',72.6,300); sb=Bosify(symmetrize(sab));
      %          sbb=deBosify(sb); subplot([ sab sbb ]);
      %
      % See also: iData_Sab/Bosify, iData_Sab/symmetrize, iData_Sqw2D
      s = iData_Sab(deBosify(iData_Sqw2D(s0), varargin{:}));
      setalias(s,'classical',   1, '[0=from measurement, with Bose factor included, 1=from MD, symmetric]');
    end
    
    function s = dynamic_range(s, varargin)
      % dynamic_range(s,Ei): crop the S(alpha,beta) to the available dynamic range
      %   for given incident neutron energy Ei [meV].
      s = iData_Sab(dynamic_range(iData_Sqw(s), varargin{:}));
      if nargout == 0
        fig=figure; h  =plot(s); set(fig, 'NextPlot','new');
      end
    end
    
    function sigma = scattering_cross_section(s, varargin)
      % scattering_cross_section(s,Ei): compute the total neutron scattering cross section for
      %   incoming neutron energy Ei [meV].
      sigma = scattering_cross_section(iData_Sqw2D(s), varargin{:});
      if nargout == 0
        fig=figure; h  =plot(sigma); set(fig, 'NextPlot','new');
      end
    end
    
    function sq = structure_factor(s)
      % structure_factor: compute the structure factor
      sq = structure_factor(iData_Sqw2D(s));
      if nargout == 0
        fig=figure; h  =plot(sq); set(fig, 'NextPlot','new');
      end
    end
    
    function [inc, multi] = incoherent(s, varargin)
      % iData_Sqw2D: incoherent: incoherent neutron scattering law estimate in the incoherent gaussian approximation
      inc = multi_phonons_incoherent(dos(s), varargin{:});
      if numel(inc) > 2, inc = inc(1:2); end
      multi = plus(inc(2:end));
      inc   = plus(inc);
      if nargout == 0
        fig=figure; 
        h  =subplot([inc multi]); 
        set(fig, 'NextPlot','new');
      end
    end
    
    function [G,multi] = multi_phonons(s, varargin)
      % iData_Sqw2D: multi_phonons: compute the integrated multi-phonon intensity from an initial density of states
      %
      % output:
      %   G:      total gDOS [p=1...]
      %   multi:  multi-phonon gDOS contribution [p=2...]
      G     = multi_phonons_dos(dos(s), varargin{:});
      multi = plus(G(2:end));
      G     = plus(G);
      if nargout == 0
        fig=figure; 
        h  =plot([ G multi ]); 
        set(fig, 'NextPlot','new');
      end
    end
    
  end % methods
  
end
