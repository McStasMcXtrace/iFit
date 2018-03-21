classdef iData_Sqw2D < iData
  % iData_Sqw2D: create a 2D S(q,w) data set (iData flavour)
  %
  % The iData_Sqw2D class is a 2D data set holding a S(q,w) dynamic structure factor
  %   aka scattering function/law.
  %   The first  axis (rows)    is the Momentum transfer [Angs-1] (wavevector). 
  %   The second axis (columns) is the Energy   transfer [meV].
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
  %
  % d   = dos(s)
  %   Compute the generalized vibrational density of states (gDOS).
  %
  % t   = thermochemistry(s)
  %   Compute and display thermochemistry quantities from the gDOS.
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
  % [inc, multi] = incoherent(s, q, T, sigma, m, n)
  %   Compute an estimate of the incoherent neutron scattering law in the incoherent gaussian approximation
  %
  % [gDOS,M]     = multi_phonons(ss)
  %   Compute the integrated multi-phonon DOS terms from an initial density of states
  %
  % input:
  %   can be a 2D iData or filename to generate a 2D Sqw object.
  %
  % output: an iData_Sqw2D object
  %
  % See also: iData, iData_Sab, iData_vDOS
  % (c) E.Farhi, ILL. License: EUPL.

  properties
  end
  
  methods
    % main instantiation method
    function obj = iData_Sqw2D(s)
      % iData_Sqw2D: create the iData_Sqw2D subclass
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
 
    end % iData_Sqw2D constructor
    
    % parameters (search for parameters in iData)
    function parameters = parseparams(s)
      % iData_Sqw2D: parseparams: search for physical quantities in object.
      % This search is also done when creating iData_Sqw2D objects.
      [s,parameters,fields] = Sqw_parameters(s);
      if nargout == 0 & length(inputname(1))
        assignin('caller',inputname(1),s);
      end
    end
    
    function f = iData(self)
      % iData_Sqw2D: iData: convert a iData_Sqw2D back to iData
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
    
    % sound_velocity ?
    % MSD ?
    % diffusion constant ?
    % compressibility ?
    % g(r) pdf
    
    function [inc, single, multi] = incoherent(s, varargin)
      % iData_Sqw2D: incoherent: incoherent neutron scattering law estimate in the incoherent gaussian approximation
      %
      % This implementation is in principle exact for an isotropic monoatomic material,
      % e.g. a liquid or powder.
      % This methodology is equivalent to the LEAPR module of NJOY ("phonon expansion")
      % to compute S(alpha,beta) from a vibrational density of states.
      % Arguments not given are searched in the data set.
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
      %   [inc, single, multi] = incoherent(sqw, q, T, sigma, m, n)
      %
      % input:
      %   sqw:the scattering law S(q,w) [iData_Sqw2D]
      %   q:  the momentum axis [Angs-1, vector]
      %   T:  temperature [K]
      %   sigma: neutron cross section [barns]
      %   m:  mass [g/mol]
      %   n:  number of iterations in the series expansion, e.g. 5
      %
      % output:
      %   inc:    total S(q,w), single+multi-phonons [iData_Sqw2D]
      %   single: so-called single-phonon incoherent S(q,w) [iData_Sqw2D]
      %   multi:  so-called multi-phonon incoherent S(q,w)  [iData_Sqw2D]
      %
      % Example:
      %   s   = iData_Sqw2D('SQW_coh_lGe.nc');
      %   inc = incoherent(s);
      %   subplot(inc);
      %
      % See also: iData_Sqw2D/multi_phonons_dos
      % (c) E.Farhi, ILL. License: EUPL.
      g   = dos(s);
      inc = incoherent(g, varargin{:});
      multi = plus(inc(3:end)); % multi-phonon
      single= plus(inc(1:2));   % elastic+1-phonon
      inc   = plus(inc);        % total
      if nargout == 0
        fig=figure; 
        h  =subplot(log10([inc single multi]),'view2'); 
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
    
    function spe = q2phi(self)
      % iData_Sqw2D: q2phi: convert a S(q,w) into a S(phi,w) iData
      spe = Sqw_q2phi(self);
      if nargout == 0
        fig=figure; 
        h  =plot(log10(s)); 
        set(fig, 'NextPlot','new');
      end
    end
    
    function s = Sab(self)
      %  iData_Sqw2D: Sab: convert a 2D S(q,w) into an S(alpha,beta). 
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
      s = Sqw2Sab(self);  % private
    end
    
    function f = saveas(self, varargin)
      % iData_Sqw2D: saveas: save S(q,w) into a file.
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
