classdef iFunc_Sqw4D < iFunc

  properties
  end
  
  methods
  
    function obj = iFunc_Sqw4D(varargin)
      % create the iFunc_Sqw4D subclass
      %
      % input:
      %   can be an iFunc or any set of parameters to generate a Sqw4D object
      obj = obj@iFunc;
      
      if nargin == 0
        m = sqw_cubic_monoatomic('defaults');
      elseif nargin == 1 && isa(varargin{1}, mfilename)
        % already an iFunc_Sqw4D
        m = varargin{1};
      elseif nargin == 1 && isa(varargin{1}, 'iFunc')
        % convert from iFunc
        m = varargin{1};
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
      for p = fieldnames(m)'
        obj.(p{1}) = m.(p{1});
      end
      obj.class = mfilename;
    end
  end % iFunc_Sqw4D constructor
  
  
end
