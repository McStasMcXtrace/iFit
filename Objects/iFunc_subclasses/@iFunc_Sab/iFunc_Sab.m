classdef iFunc_Sab < iFunc
  % iFunc_Sab: create an iFunc_Sab from e.g. an iFunc 2D object
  %
  % The iFunc_Sab class is a 2D model holding a S(alpha,beta) dynamic structure factor.
  % The first axis is alpha [e.g unitless moment], 2nd is beta [e.g. unitless energy].
  %
  % We commonly define (h stands for hbar):
  %     alpha=  h2q2/2MkT = (Ei+Ef-2*mu*sqrt(Ei*Ef))/AkT   unit-less moment
  %     beta = -hw/kT     = (Ef-Ei)/kT                     unit-less energy
  %     A    = M/m                
  %     mu   = cos(theta) = (Ki.^2 + Kf.^2 - q.^2) ./ (2*Ki.*Kf)
  %     m    = mass of the neutron
  %     M    = mass of the scattering target [g/mol]
  %
  % The 'quantum' dynamic structure factor S(a,b) should obey the following rules:
  %
  %   \int S(a,b) db   = 1      for a pure incoherent scatterer
  %   \int S(a,b) b db = -a
  %   S(a,b) = exp(-b) S(a,-b)  so-called detailed balance
  %
  % The symmetric/classical S*(a,b) can be defined as:
  %   S*(a,b) = exp(b/2) S(a,b) 
  % but other 'quantum corrections' are possible.
  %
  % Useful methods for this iFunc flavour:
  %
  % methods(iFunc_Sab)
  %   all iFunc methods can be used.
  % iFunc_Sab(a)
  %   convert a 2D model [a=iFunc class] into an iFunc_Sab to give access to
  %   the methods below.
  %
  % input:
  %   can be an iFunc or struct or any set of parameters to generate a Sab object.
  %
  % output: an iFunc_Sab object
  %
  % See also: sqw_phonons, sqw_cubic_monoatomic, iFunc

  properties
  end
  
  methods
  
    function obj = iFunc_Sab(varargin)
      
      obj = obj@iFunc;
      
      if nargin == 0
        return
      elseif nargin == 1 && isa(varargin{1}, mfilename)
        % already an iFunc_Sab
        m = varargin{1};
      elseif nargin == 1 && isa(varargin{1}, 'iFunc_Sqw2D')
        m = Sqw2Sab(varargin{1});
      elseif nargin == 1 && isa(varargin{1}, 'iFunc')
        % convert from iFunc
        m = varargin{1};
      elseif nargin == 1 && isa(varargin{1}, 'struct')
        m = iFunc(varargin{1});
      end
      
      % from here, we must have either an iFunc or an iFunc_Sab
      if isa(m, mfilename)
        obj = m; obj.class = mfilename;
        return;
      elseif ~isa(m, 'iFunc')
        error([ mfilename ': the given input ' class(m) ' does not seem to be convertible to iFunc_Sab.' ])
      end
      
      % check if the Sab subclass is appropriate
      flag = false;
      if ndims(m) == 2 % must be S(q,w) or S(a,b)
        flag = true;
      end
      if ~flag
        error([ mfilename ': the given iFunc model does not seem to be an Sab flavour object.' ])
      end
      
      % transfer properties
      % this is a safe way to instantiate a subclass
      warning off MATLAB:structOnObject
      m = struct(m);
      for p = fieldnames(m)'
        obj.(p{1}) = m.(p{1});
      end
      obj.class = mfilename;
    end % iFunc_Sab constructor
    
    function f = iFunc(self)
      % iFunc_Sab: convert a single iFunc_Sab back to iFunc
      f = iFunc;
      warning off MATLAB:structOnObject
      self = struct(self);
      for p = fieldnames(self)'
        f.(p{1}) = self.(p{1});
      end
      f.class = 'iFunc';
    end
    
    function f = iData(self, varargin)
      % iFunc_Sab: iData: evaluate a 2D Sab into an iData object
      %
      %   iData(self, p, q, w)
      
      % check for alpha beta grid
      if isempty(varargin),  varargin{end+1} = []; end % parameters
      if numel(varargin) <2, varargin{end+1} = linspace(0,150,31); end
      if numel(varargin) <3, varargin{end+1} = linspace(0.1,400,41); end
      s = iFunc(self);
      f = transpose(iData(s,varargin{:}));
      xlabel(f, 'Alpha [h^2q^2/2MkT]');
      ylabel(f, 'Beta [-hw/kT]');
      title(f, self.Name);
      self.UserData = s.UserData;
      if ~isempty(inputname(1))
        assignin('caller',inputname(1),self); % update in original object
      end
    end
    
    function [signal, self, ax, name] = feval(self, varargin)
      % iFunc_Sab: feval: evaluate the Model on alpha/beta grid
      if isempty(varargin),  varargin{end+1} = []; end % parameters
      if numel(varargin) <2, varargin{end+1} = linspace(0,150,31); end
      if numel(varargin) <3, varargin{end+1} = linspace(0.1,400,41); end
      [signal, self, ax, name] = feval@iFunc(self, varargin{:});
      if ~isempty(inputname(1))
        assignin('caller',inputname(1),self); % update in original object
      end
    end
    
    function h = plot(self, varargin)
      % iFunc_Sab: plot: plot the S(alpha,beta) model
      h = plot@iFunc(self, varargin{:});
      ylabel('Alpha=h2q2/2MkT [1]');
      xlabel('Beta=-hw/kT [1]');
    end
    
    
  end % methods
  
end % classdef
  

% private functions used in the class ----------------------------------------


