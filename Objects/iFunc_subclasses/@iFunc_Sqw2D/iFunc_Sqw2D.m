classdef iFunc_Sqw2D < iFunc
  % iFunc_Sqw2D: create an iFunc_Sqw2D from e.g. an iFunc 2D object
  %
  % The iFunc_Sqw2D class is a 2D model holding a S(q,w) dynamic structure factor.
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
  % b = band_structure(a)
  %   Compute the band structure
  % d = dos(a)
  %   Compute the density of states
  % m = max(a)
  %   Evaluate the maximum S(q,w) dispersion energy
  % publish(a)
  %   Generate a readable document with all results
  % t = thermochemistry(a)
  %   Compute and display thermochemistry quantities
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
        m = sqw_cubic_monoatomic('defaults');
      elseif nargin == 1 && isa(varargin{1}, mfilename)
        % already an iFunc_Sqw2D
        m = varargin{1};
      elseif nargin == 1 && isa(varargin{1}, 'iFunc')
        % convert from iFunc
        m = varargin{1};
      elseif nargin == 1 && isa(varargin{1}, 'struct')
        m = iFunc(varargin{1});
      else
        % create new Sqw2D model and convert it to iFunc_Sqw2D
        m = sqw_phonons(varargin{:}); 
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
      if numel(varargin) <2, varargin{end+1} = linspace(0,0.5,30); end
      if numel(varargin) <3, varargin{end+1} = linspace(0.01,max(self)*1.2,11); end
      s = iFunc(self);
      f =iData(s,varargin{:});
      xlabel(f, 'Q [Angs]');
      ylabel(f, 'Energy [meV]');
      title(f, self.Name);
      self.UserData = s.UserData;
      if ~isempty(inputname(1))
        assignin('caller',inputname(1),self); % update in original object
      end
    end
    
    function [fig, s, k] = plot(self, varargin)
      % iFunc_Sqw2D: plot: plot dispersions along principal axes and gDOS
      if isempty(varargin), varargin = { 'plot meV' }; end
      [s,k,fig]=band_structure(self, varargin{:});
      if ~isempty(inputname(1))
        assignin('caller',inputname(1),self); % update in original object
      end
    end % plot
  
    % methods for Sqw 2D
    
    % we evaluate the model, extend its Expression:
    % convert it to an iData_Sqw2D
    % apply one of the following method:
    
    % Bosify
    % deBosify
    % dos/gdos -> 1D iFunc
    % structure_factor/sq -> 1D iFunc
    % thermochemistry -> 1D iFunc array
    % moments -> 1D iFunc array
    % dynamic_range
    % scattering_cross_section -> 1D iFunc
    % incoherent (->iData_Sqw2D->incoherent)
    % coherent (+sq)
    % multi_phonons
    % symmetrize
    % qw2phiw
    % qw2qt
    % qw2phi
    % Sqw2ddcs
    % ddcs2Sqw
    % Sab
    
    
  end % methods
  
end % classdef
  

% private functions used in the class ----------------------------------------


