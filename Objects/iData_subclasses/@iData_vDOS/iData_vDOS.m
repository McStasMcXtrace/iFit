classdef iData_vDOS < iData

  % iData_vDOS: create a vibrational density of states data set (iData flavour)
  %
  % The iData_vDOS class is a 1D data set holding a vibrational density of states
  % often labelled as g(w)
  %
  % Example: g=iData_vDOS('SQW_coh_lGe.nc')
  %
  % Useful methods for this iData flavour:
  %
  % methods(iData_Sqw2D)
  %   all iData methods can be used.
  % iData_vDOS(g)
  %   convert input [e.g. a 2D or 4D Sqw, or a 1D vDOS] into an iData_vDOS to give access to
  %   the methods below.
  %
  % d = dos(g)
  %   Compute/get the vibrational density of states (vDOS).
  %
  % t = thermochemistry(g)
  %   Compute and display thermochemistry quantities from the vDOS.
  %
  % [Sqw,Iqt,DW] = multi_phonons_incoherent(g, q, T, sigma, m, n)
  %   Compute the so-called multi-phonon dynamic structure factor expansion in 
  %   the incoherent gaussian approximation. The neutron scattering incoherent 
  %   gaussian dynamic structure factor is computed, as well as the corresponding
  %   intermediate scattering function, and the Debye-Waller factor.
  %
  % [Gw, Tsym] = multi_phonons_vDOS(g, Ki, T, sigma, m, phi, n)
  %   Compute the so-called multi-phonon density of states expansion in the 
  %   neutron scattering incoherent gaussian approximation.
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
    function obj = iData_vDOS(s)
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
    
    function d = dos(s)
    % iData_vDOS: dos: return the object. 
    %   For easier syntax and compatibility with other classes.
      d = s;
    end
    
    % parameters (search for parameters in iData)
    function parameters = parseparams(s)
      % iData_vDOS: parseparams: search for physical quantities in object.
      % This search is also done when creating iData_vDOS objects.
      [s,parameters,fields] = Sqw_parameters(s);
    end
    
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
    end
    
  end
  
end
