classdef iData_vDOS < iData

  % iData_vDOS: create a vibrational density of states data set (iData flavour)
  %
  % The iData_vDOS class is a 1D data set holding a density of states
  %   often labelled as g(w). This usually can store a 'true' vibrational density 
  %   of states (vDOS), or a generalised density of states (gDOS).
  %
  % Example: g=iData_vDOS('SQW_coh_lGe.nc')
  %
  % Useful methods for this iData flavour:
  %
  % methods(iData_vDOS)
  %   all iData methods can be used.
  % iData_vDOS(g)
  %   convert input [e.g. a 2D or 4D Sqw, or a 1D vDOS] into an iData_vDOS to give access to
  %   the methods below.
  %
  % d = dos(g)
  %   Get the vibrational density of states (vDOS).
  %
  % t = thermochemistry(g)
  %   Compute and display thermochemistry quantities from the vDOS.
  %
  % [Sqw,Iqt,DW] = incoherent(g, q, T, sigma, m, n, DW)
  %   Compute the incoherent gaussian approximation scattering law S(q,w) from an initial 
  %   density of states vDOS.
  %
  % [Gw, Tsym] = multi_phonons(g, Ki, T, sigma, m, phi, n)
  %   compute the integrated multi-phonon generalised density of states (gDOS) 
  %   from an initial vibrational density of states (vDOS) in the neutron 
  %   scattering incoherent gaussian approximation.
  %   The generalised density of states (gDOS) is the sum of the terms in this expansion.
  %
  % gDOS = gdos(vDOS, Ki, T, sigma, m, phi, n)
  %   compute the generalised density of states (gDOS) from a vibrational 
  %   density of states (vDOS) in the neutron scattering incoherent gaussian approximation.
  %
  % vDOS = vdos(gDOS, Ki, T, sigma, m, phi, n)
  %   compute the 'true' vibrational density of states from a gDOS estimate
  %   in the neutron scattering incoherent gaussian approximation.
  %
  % input:
  %   can be an iData (1D vDOS, 2D Sqw, 4D Sqw) or filename to generate a 1D vDOS object.
  %
  % output: an iData_vDOS object
  %
  % See also: iData

  properties
  end
  
  methods
  
    % main instantiation method
    function obj = iData_vDOS(s, varargin)
      % iData_vDOS: create the iData_vDOS subclass
      %
      % input:
      %   must be any of
      %    1D iData object [energy,density of states] 
      %    2D iData S(q,w), S(phi,t), S(alpha,beta) or S(phi,w) [iData_Sqw2D]
      %    4D iData S(q,w) [iData_Sqw4D]
      %    filename containing any of the above
      %
      % Example: s=iData_vDOS('SQW_coh_lGe.nc');
      
      % could also support
      %    2D iFunc model holding a S(q,w)
      %    4D iFunc model holding a S(q,w)
      
      obj = obj@iData;
      obj.class = mfilename;
      if ~nargin, return; end
      if ~isa(s, 'iData')
        s = iData(s, varargin{:});
      end
      
      % convert/test
      m = [];
      if isa(s, 'iData')
        if ndims(s) == 1
          m = s;
        elseif ndims(s) == 2
          m = dos(iData_Sqw2D(s));
        elseif ndims(s) == 4
          m = dos(powder(iData_Sqw4D(s)));
        end
      end
      
      if isempty(m)
        error([ mfilename ': the given input ' class(s) ' does not seem to be convertible to iData_vDOS.' ])
      end
      
      m = Sqw_parameters(m);
      
      % transfer properties
      % this is a safe way to instantiate a subclass
      warning off MATLAB:structOnObject
      m = struct(m);
      for p = fieldnames(m)'
        obj.(p{1}) = m.(p{1});
      end
      obj.class = mfilename;
 
    end % iData_vDOS constructor
    
    function d = dos(self)
    % iData_vDOS: dos: return the object. 
    %   For easier syntax and compatibility with other classes.
      d = self;
    end % dos
    
    function g = gdos(self, varargin)
    % iData_vDOS: gdos: compute the generalised density of states from a vibrational density of states. Should be compared to experimental gDOS.
      g   = multi_phonons(self, varargin{:});
      % normalise the 1-phonon term to 1 and apply to others
      nrm = trapz(g(1));
      g   = g ./ nrm; % normalise to gDOS 1-phonon term
      g   = plus(g); % effective gDOS == \int q S(q,w) dq
    end % gdos
    
    % parameters (search for parameters in iData)
    function parameters = parseparams(s)
      % iData_vDOS: parseparams: search for physical quantities in object.
      % This search is also done when creating iData_vDOS objects.
      [s,parameters,fields] = Sqw_parameters(s);
    end % parseparams
    
    function f = iData(self)
      % iData_vDOS: iData: convert a iData_vDOS back to iData
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
    end % iData
    
    function v = vdos(self, varargin)
      % iData_vDOS: vdos: compute the 'true' vibrational density of states from a gDOS estimate.
      %
      % Example:
      %   gdos = dos(iData_Sqw2D)
      %   vdos = vdos(g); % refine and remove multi-phonons
      %
      % input:
      %   self:   a generalised density of states (gDOS), e.g. \int q S(q,w) dq
      %   order:  restrict the multi-phonon expansion to given order
      %
      % output:
      %   v:      vDOS
      
      % self is e.g. \int q S(q,w) dq = \int S(theta,w) sin(theta) d(theta)
      
      % we want that self == gdos(v) == sum(multi_phonons terms)
      
      % norm(gdos) > 1. norm(v) == 1
      
      self = self./trapz(self); % normalise initial gDOS to 1. In fact norm should be higher if contains multi-phonons
      
      v = self;         % we start by assuming the gDOS is not too far from the vdos
      
      flag = false;
      
      while ~flag
        % compute the gDOS from the vDOS
        g   = gdos(v, varargin{:});
        nrm = trapz(g)  % total norm of the gDOS, including multi-phonons, nrm > 1
        
        % we compare the initial gDOS (self) to our current estimate (g)
        % our initial gDOS has norm 1, we want it to have norm 'nrm' for comparison
        ratio = self*nrm./g;
        
        % check if we have converged (correction is small) e.g. self == gdos(v)
        flag = all(abs(ratio - 1) < .01); % changes smaller than 1 percent everywhere
        
        % correct our vDOS estimate. When ratio > 1, self is too large, and we should increase 'v'
        v = v .* ratio;
      end
      
    end % vdos
    
    
  end % methods
  
end % class
