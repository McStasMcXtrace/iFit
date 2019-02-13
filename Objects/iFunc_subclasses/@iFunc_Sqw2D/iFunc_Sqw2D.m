classdef iFunc_Sqw2D < iFunc
  % iFunc_Sqw2D: create an iFunc_Sqw2D from e.g. an iFunc 2D object
  %
  % The iFunc_Sqw2D class is a 2D model holding a S(q,w) dynamic structure factor.
  % The first axis is Q [e.g Angs-1], 2nd is Energy [e.g. meV].
  %
  % The 'quantum' dynamic structure factor S(q,w) should obey the following rules:
  %   S(q) = \int S(q,w) dw = <S(q,w)>                 structure factor [moment 0]
  %   Er   = \int w*S(q,w) dw = <wS(q,w)> = h2q2/2M       recoil energy [moment 1]
  %   S(q,-w) = exp(-hw/kT) S(q,w)                      so-called detailed balance
  %
  % The symmetric/classical S*(q,w) can be defined as:
  %   S*(q,w) = exp(hw/2kT) S(q,w)
  % but other 'quantum corrections' are possible.
  %
  % Useful methods for this iFunc flavour:
  %
  % methods(iFunc_Sqw2D)
  %   all iFunc methods can be used.
  % iFunc_Sqw2D(a)
  %   convert a 2D model [a=iFunc class] into an iFunc_Sqw2D to give access to
  %   the methods below.
  % plot(a)
  %   plot the band structure and density of states
  %
  % input:
  %   can be an iFunc or struct or any set of parameters to generate a Sqw2D object.
  %   when not given an iFunc, the parameters to sqw_phonons are expected.
  %
  % output: an iFunc_Sqw2D object
  %
  % See also: sqw_phonons, sqw_cubic_monoatomic, iFunc

  properties
  end
  
  methods
  
    function obj = iFunc_Sqw2D(varargin)
      
      obj = obj@iFunc;
      
      if nargin == 0
        return
      elseif nargin == 1 && isa(varargin{1}, mfilename)
        % already an iFunc_Sqw2D
        m = iFunc(varargin{1}); % check
      elseif nargin == 1 && isa(varargin{1}, 'iFunc')
        % convert from iFunc
        m = varargin{1};
      elseif nargin == 1 && isa(varargin{1}, 'struct')
        m = iFunc(varargin{1});
      end
      
      % from here, we must have either an iFunc or an iFunc_Sqw2D
      if isa(m, mfilename)
        obj = m; obj.class = mfilename;
        return;
      elseif ~isa(m, 'iFunc')
        error([ mfilename ': the given input ' class(m) ' does not seem to be convertible to iFunc_Sqw2D.' ])
      end
      
      % check if the Sqw2D subclass is appropriate
      flag = false;
      if ndims(m) == 2 % must be S(hkl,w)
        flag = true;
      end
      if ~flag
        error([ mfilename ': the given iFunc model does not seem to be an Sqw2D flavour object.' ])
      end
      
      % transfer properties
      % this is a safe way to instantiate a subclass
      warning off MATLAB:structOnObject
      m = struct(m);
      for p = fieldnames(m)'
        obj.(p{1}) = m.(p{1});
      end
      obj.class = mfilename;
    end % iFunc_Sqw2D constructor
    
    function f = iFunc(self)
      % iFunc_Sqw2D: convert a single iFunc_Sqw2D back to iFunc
      f = iFunc;
      warning off MATLAB:structOnObject
      self = struct(self);
      for p = fieldnames(self)'
        f.(p{1}) = self.(p{1});
      end
      f.class = 'iFunc';
    end
    
    function f = iData(self, varargin)
      % iFunc_Sqw2D: iData: evaluate a 2D Sqw into an iData object
      %
      %   iData(self, p, q, w)
      
      % check for Q W grid
      if isempty(varargin),  varargin{end+1} = []; end % parameters
      if numel(varargin) <2, varargin{end+1} = linspace(0,3,31); end
      if numel(varargin) <3, varargin{end+1} = linspace(-20,20,41); end
      s = iFunc(self);
      f = transpose(iData(s,varargin{:}));
      xlabel(f, 'Q [Angs]');
      ylabel(f, 'Energy [meV]');
      title(f, self.Name);
      self.UserData = s.UserData;
      if ~isempty(inputname(1))
        assignin('caller',inputname(1),self); % update in original object
      end
    end
    
    function [signal, self, ax, name] = feval(self, varargin)
      % iFunc_Sqw2D: feval: evaluate the Model on QW grid
      if isempty(varargin),  varargin{end+1} = []; end % parameters
      if numel(varargin) <2, varargin{end+1} = linspace(0,3,31); end
      if numel(varargin) <3, varargin{end+1} = linspace(-20,20,41); end
      [signal, self, ax, name] = feval@iFunc(self, varargin{:});
      if ~isempty(inputname(1))
        assignin('caller',inputname(1),self); % update in original object
      end
    end
    
    function h = plot(self, varargin)
      % iFunc_Sqw2D: plot: plot the S(q,w) model
      h = plot@iFunc(self, varargin{:});
      ylabel('Q [Angs-1]');
      xlabel('Energy [meV]');
    end
    
    function spw = angle(self)
      % iFunc_Sqw2D: angle: convert a S(q,w) into a S(phi,w) Model (scattering angle)
      %
      % convert: iFunc_Sqw2D S(q,w) -> S(phi,w)
      %
      % New model parameters:
      %       Ei     Incident neutron energy [meV]
      %
      %  The incident neutron energy can be computed using:
      %   Ei = 2.0721*Ki^2 = 81.8042/lambda^2 with Ki in [Angs-1] and lambda in [Angs]
      %
      % spw = angle(s)
      %
      % input:
      %   s:      Sqw 2D model with q as 1st axis (Angs-1), w as 2nd axis (meV).
      %
      % output:
      %   spw:    S(phi,w) 2D model [iFunc, phi in deg]
      spw = qw2phiw(self);
    end
    
    function sab = Sab(self)
      % iFunc_Sqw2D: Sab: convert a S(q,w) into a S(alpha,beta) Model
      %
      % convert: iFunc_Sqw2D S(q,w) -> S(alpha,beta)
      %
      % New model parameters:
      %       Temperature                      [K]
      %       Mass       Material molar weight [g/mol]
      
      sab = iFunc_Sab(self);
    end % Sab
  
    % methods for iFunc Sqw 2D, similar to the iData Sqw2D ones
    
    % methods that retain the type (q,w)
    %   incoherent(q, T, m, n, DW, vDOS ) using dos.incoherent, hard work. Use iData_vDOS ?
    
    % methods that change the type (q,w) -> something else (iFunc 1D)
    %   multi_phonons -> 1D(w)
    %   scattering_cross_section(Ei_min, Ei_max, Ei_n, Mass) -> 1D iFunc(Ei) -> new axis !!
    %   moments(M,T, classical) -> 1D iFunc array ? from dos
    %   structure_factor/sq -> 1D iFunc
    %   thermochemistry(T) -> 1D iFunc array ? from dos
    
    % methods that change the type (q,w) -> something else (iFunc 2D)
    %   Sab(M,T)
    
  end % methods
  
end % classdef
  

% private functions used in the class ----------------------------------------


