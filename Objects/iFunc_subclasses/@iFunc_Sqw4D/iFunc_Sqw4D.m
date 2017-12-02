classdef iFunc_Sqw4D < iFunc
  % iFunc_Sqw4D: create an iFunc_Sqw4D from e.g. an iFunc 4D object
  %
  % The iFunc_Sqw4D class is a 4D model holding a S(q,w) dynamic structure factor.
  %
  % Example: s=sqw_cubic_monoatomic
  %
  % Useful methods for this iFunc flavour:
  %
  % methods(iFunc_Sqw4D)
  %   all iFunc methods can be used.
  % iFunc_Sqw4D(a)
  %   convert a 4D model [a=iFunc class] into an iFunc_Sqw4D to give access to
  %   the methods below.
  % plot(a)
  %   plot the band structure and density of states
  % plot3(a)
  %   plot the S(q,w) dispersion in the QH=0 plane in 3D
  % scatter3(a)
  %   plot the S(q,w) dispersion in the QH=0 plane in 3D as coloured points
  % slice(a)
  %   Slice and isosurface volume exploration of the S(q,w) dispersion in the QH=0 plane
  % b = band_structure(a)
  %   Compute the band structure
  % d = dos(a)
  %   Compute the density of states
  % m = max(a)
  %   Evaluate the maximum S(q,w) dispersion energy
  % powder(a)
  %   Compute the powder average of the 4D S(q,w) dispersion
  % publish(a)
  %   Generate a readable document with all results
  % t = thermochemistry(a)
  %   Compute and display thermochemistry quantities
  %
  % input:
  %   can be an iFunc or struct or any set of parameters to generate a Sqw4D object.
  %   when not given an iFunc, the parameters to sqw_phonons are expected.
  %
  % output: an iFunc_Sqw4D object
  %
  % See also: sqw_phonons, sqw_cubic_monoatomic, iFunc

  properties
  end
  
  methods
  
    function obj = iFunc_Sqw4D(varargin)
      
      obj = obj@iFunc;
      
      if nargin == 0
        m = sqw_cubic_monoatomic('defaults');
      elseif nargin == 1 && isa(varargin{1}, mfilename)
        % already an iFunc_Sqw4D
        m = varargin{1};
      elseif nargin == 1 && isa(varargin{1}, 'iFunc')
        % convert from iFunc
        m = varargin{1};
      elseif nargin == 1 && isa(varargin{1}, 'struct')
        m = iFunc(varargin{1});
      else
        % create new Sqw4D model and convert it to iFunc_Sqw4D
        m = sqw_phonons(varargin{:}); 
      end
      
      % from here, we must have either an iFunc or an iFunc_Sqw4D
      if isa(m, mfilename)
        obj = m; obj.class = mfilename;
        return;
      elseif ~isa(m, 'iFunc')
        error([ mfilename ': the given input ' class(m) ' does not seem to be convertible to iFunc_Sqw4D.' ])
      end
      
      % check if the Sqw4D subclass is appropriate
      flag = false;
      if ndims(m) == 4 % must be S(hkl,w)
        flag = true;
      end
      if ~flag
        error([ mfilename ': the given iFunc model does not seem to be an Sqw4D flavour object.' ])
      end
      
      % transfer properties
      % this is a safe way to instantiate a subclass
      warning off MATLAB:structOnObject
      m = struct(m);
      for p = fieldnames(m)'
        obj.(p{1}) = m.(p{1});
      end
      obj.class = mfilename;
    end % iFunc_Sqw4D constructor
    
    function f = iFunc(self)
      % iFunc_Sqw4D: convert a single iFunc_Sqw4D back to iFunc
      f = iFunc;
      self = struct(self);
      for p = fieldnames(self)'
        f.(p{1}) = self.(p{1});
      end
      f.class = 'iFunc';
    end
    
    function f = iData(self, varargin)
      % iFunc_Sqw4D: iData: evaluate a 4D Sqw into an iData object
      f =iData(iFunc(self),varargin{:});
      xlabel(f, 'QH [rlu]');
      ylabel(f, 'QK [rlu]');
      zlabel(f, 'QL [rlu]');
      clabel(f, 'Energy [meV]');
      title(f, self.Name);
    end
    
    function [fig, s, k] = plot(self, varargin)
      % iFunc_Sqw4D: plot: plot dispersions along principal axes and vDOS
      if isempty(varargin), varargin = { 'plot meV' }; end
      [s,k,fig]=band_structure(self, varargin{:});
      if ~isempty(inputname(1))
        assignin('caller',inputname(1),self); % update in original object
      end
    end % plot
    
    function [h, f] = plot3(s, varargin)
      % iFunc_Sqw4D: plot3: plot a 3D view of the dispersions in H=0 plane
      f = feval_fast(s);
      % plot in 3D
      h = plot3(f, varargin{:}); % h=0
      if ~isempty(inputname(1))
        assignin('caller',inputname(1),s); % update in original object
      end
    end % plot3
    
    function [h, f] = scatter3(s, varargin)
      % iFunc_Sqw4D: scatter3: plot a 3D scatter view of the dispersions in H=0 plane
      f = feval_fast(s);
      % plot in 3D
      h = scatter3(f, varargin{:}); % h=0
      if ~isempty(inputname(1))
        assignin('caller',inputname(1),s); % update in original object
      end
    end % scatter3
    
    function [h, f] = slice(s, varargin)
      % iFunc_Sqw4D: slice: plot an editable 3D view of the dispersions in H=0 plane
      f = feval_fast(s);
      % plot in 3D
      h = slice(f, varargin{:}); % h=0
      if ~isempty(inputname(1))
        assignin('caller',inputname(1),s); % update in original object
      end
    end % slice
  
    % methods for Sqw 4D
    
    % bosify
    % debosify
    % gdos
    % sq
  end % methods
  
end % classdef
  

% private functions used in the class ----------------------------------------


