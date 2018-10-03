classdef iFunc_Sqw2D < iFunc
  % iFunc_Sqw2D: create an iFunc_Sqw2D from e.g. an iFunc 2D object
  %
  % The iFunc_Sqw2D class is a 2D model holding a S(q,w) dynamic structure factor.
  % The first axis is Q [e.g Angs-1], 2nd is Energy [e.g. meV].
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
        return
      elseif nargin == 1 && isa(varargin{1}, mfilename)
        % already an iFunc_Sqw2D
        m = varargin{1};
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
      if numel(varargin) <3, varargin{end+1} = linspace(0.01,20,21); end
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
  
    % methods for Sqw 2D
    
    % we evaluate the model, extend its Expression.
    %   convert it to an iData_Sqw2D and apply one of the following method:
    
    % methods that retain the type (q,w)
    %   Bosify(T, method)
    %   deBosify(T, method)
    %   dynamic_range(Ei, angle_min, angle_max)
    %   incoherent(q, T, m, n, DW)
    %   coherent (iData sq or d-spacing value)
    %   symmetrize
    %   Sqw2ddcs(Ei)
    %   ddcs2Sqw(Ei)
    
    % when changing type, the new object must compute back the (q,w) axes, then 
    % evaluate the original iFunc_Sqw2D model, and apply the final conversion with
    % 'new' axes.
    
    % methods that change the type (q,w) -> something else (iFunc 1D)
    %   dos/gdos(T, DW, method) -> 1D iFunc(w)
    %   multi_phonons -> 1D(w)
    %   scattering_cross_section(Ei_min, Ei_max, Ei_n, Mass) -> 1D iFunc(Ei) -> new axis !!
    %   moments(M,T, classical) -> 1D iFunc array ?
    %   structure_factor/sq -> 1D iFunc
    %   thermochemistry(T) -> 1D iFunc array ?
    %   muphocor(T, amasi, sigi, conci) -> 1D array
    
    % methods that change the type (q,w) -> something else (iFunc 2D)
    %   Sab(M,T)
    %   qw2phiw(lambda)
    %   qw2qt(lambda, chwidth, dist)
    %   qw2phi(lambda)
    
    function f1 = addtoexpr(f0, method, varargin)
      % addtoexpr: catenate a given iData_Sqw2D method to the iFunc object Expression
      
      % all parameters must be given explicitly, NOT from internal physical parameter search
      % varargin = {'pars,value, ...}
      %   add new parameters when value is scalar
      %   add to UserData when not scalar (and get it back for eval of iData_Sqw2D method)
      
      % varargin must be stored in the UserData ? or added as parameters for single values ?
      % this depends on the iData_Sqw2D method
      % what about parseparams and Sqw_check for parameters ? -> iFunc parameters OK
      
      % we use the '+' with a char string to catenate the expression
      f1            = f0 + [ ...
        'q=x; w=y; signal = iData_Sqw2D(iData(q,w,signal));' ...
        'signal = ' method '(signal, varargin);' ...
        'signal = getaxis(signal, 0);' ...
        ];

    end
    
    
  end % methods
  
end % classdef
  

% private functions used in the class ----------------------------------------


