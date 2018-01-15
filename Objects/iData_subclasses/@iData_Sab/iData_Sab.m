classdef iData_Sab < iData
  % iData_Sab: create a 2D S(alpha,beta) data set (iData flavour)
  %
  % The iData_Sab class is a 2D data set holding a S(alpha,beta) dynamic structure factor.
  %
  %  The S(alpha,beta) is a representation of the dynamic structure factor 
  %  using unitless momentum and energy variables defined as:
  %     alpha= h2q2/2MkT  = (Ei+Ef-2*mu*sqrt(Ei*Ef))/AkT
  %     beta = -hw/kT     = (Ef-Ei)/kT
  %     A    = M/m
  %     mu   = cos(theta) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf)
  %
  % Example: s=iData_Sab('SQW_coh_lGe.nc')
  %
  % Useful methods for this iData flavour:
  %
  % methods(iData_Sab)
  %   all iData methods can be used.
  % iData_Sab(s)
  %   convert input [e.g. a 2D iData object] into an iData_Sab to give access to
  %   the methods below.
  %
  % d   = dos(s)
  %   Compute the generalized vibrational density of states (gDOS).
  %
  % t   = thermochemistry(s)
  %   Compute and display thermochemistry quantities from the gDOS.
  %
  % m   = moments(s)
  %   Compute the S(alpha,beta) moments/sum rules (harmonic frequencies).
  %
  % sym = symmetrize(s)
  %   Extend the S(alpha,beta) in both 'beta' sides.
  %
  % sb  = Bosify(s)
  %   Apply the 'Bose' factor (detailed balance) to a classical data set.
  %
  % s   = deBosify(sb)
  %   Remove Bose factor (detailed balance) from an 'experimental/quantum' data set.
  %
  % p   = parseparams(s)
  %   Search for physical quantities in S(alpha,beta) data set.
  %
  % sqw = Sqw(sab)
  %   Convert an S(alpha,beta) to an S(q,w), which is roughly independent of the temperature.
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
  % input:
  %   can be an iData or filename to generate a 2D Sab object.
  %
  % output: an iData_Sab object
  %
  % See also: iData
  
  properties
  end
  
  methods
    % main instantiation method
    function obj = iData_Sab(s, varargin)
      % iData_Sab: create the iData_Sab subclass
      %
      % input:
      %   must be a 2D iData object holding a S(q,w), S(phi,t), S(alpha,beta) or S(phi,w)
      %   or file name.
      %
      % Example: s=iData_Sab('SQW_coh_lGe.nc');
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
      % iData_Sab: parseparams: search for physical quantities in object.
      % This search is also done when creating iData_Sab objects.
      [s,parameters,fields] = Sqw_parameters(s);
    end
    
    function f = iData(self)
      % iData_Sab: iData: convert a iData_Sab back to iData
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
    
    function [DOS, fig] = dos(self)
      % iData_Sab: dos: compute the generalised density of states (gDOS)
      DOS = dos(iData_Sqw2D(self));
      
      if nargout == 0 && ~isempty(DOS)
        fig=figure; 
        % plot total DOS
        h=plot(DOS); set(h,'LineWidth',2);
        set(fig, 'NextPlot','new');
      else fig = [];
      end
    end
    
    function t = thermochemistry(s, T, options)
      % iData_Sab: thermochemistry: compute thermodynamic quantities for 2D S(alpha,beta) data sets.
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
    
    function sigma=moments(data, varargin)
      % iData_Sab: moments=moments(sqw, M, T, classical): compute Sqw moments/sum rules (harmonic frequencies)
      sigma = moments(iData_Sqw2D(data), varargin{:});
    end
    
    function sym = symmetrize(s)
      % iData_Sab: symmetrize(s): extend the S(alpha,beta) in both alpha/energy sides
      s = symmetrize(iData_Sqw2D(s));
    end
    
    function s = Bosify(s0, varargin)
      % Bosify: apply the 'Bose' factor (detailed balance) to a classical data set.
      s = iData_Sab(Bosify(iData_Sqw2D(s0), varagin{:}));
      setalias(s,'classical',   0, '[0=from measurement, with Bose factor included, 1=from MD, symmetric]');
    end
    
    function s = deBosify(s0, varargin)
      % deBosify: remove Bose factor (detailed balance) from an 'experimental' data set.
      s = iData_Sab(deBosify(iData_Sqw2D(s0), varagin{:}));
      setalias(s,'classical',   1, '[0=from measurement, with Bose factor included, 1=from MD, symmetric]');
    end
    
    function s = dynamic_range(s, varargin)
      % dynamic_range(s,Ei): crop the S(alpha,beta) to the available dynamic range
      %   for given incident neutron energy.
      s = iData_Sab(dynamic_range(iData_Sqw(s), varargin{:}));
    end
    
    function sigma = scattering_cross_section(s, varargin)
      % scattering_cross_section(s,Ei): compute the total neutron scattering cross section for
      %   incoming neutron energy.
      sigma = scattering_cross_section(iData_Sqw2D(s), varargin{:});
    end
    
    function s = structure_factor(s)
      % structure_factor: compute the structure factor
      s = structure_factor(iData_Sqw2D(s));
    end
  end % methods
  
end
