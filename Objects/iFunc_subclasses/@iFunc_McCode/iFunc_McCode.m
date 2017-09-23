classdef iFunc_McCode < iFunc

  properties
  end

  methods
    function obj = iFunc_McCode(varargin)
      % create the subclass
      obj = obj@iFunc;
      if nargin == 0
        return
      elseif nargin == 1 && isa(varargin{1}, mfilename)
        % already an iFunc_McCode
        m = varargin{1};
      elseif nargin == 1 && isa(varargin{1}, 'iFunc')
        % convert from iFunc
        m = varargin{1};
      else
        % create new McCode model and convert it to iFunc_McCode
        m = mccode(varargin{:}); 
      end
      
      % from here, we must have either an iFunc or an iFunc_McCode
      if isa(m, mfilename)
        obj = m; return;
      elseif ~isa(m, 'iFunc')
        error([ mfilename ': the given input ' class(m) ' does not seem to be convertible to iFunc_McCode.' ])
      end
      
      % check if the McCode subclass is appropriate
      flag = false;
      for t = { 'McCode','McStas',' McXtrace' }'
        if strfind(lower(m.Name),lower(t{1}))
          flag = true; break;
        end
      end
      if ~flag
        error([ mfilename ': the given iFunc model does not seem to be an McCode flavour object.' ])
      end
      
      % transfer properties
      % this is a safe way to instantiate a subclass
      for p = fieldnames(m)'
        obj.(p{1}) = m.(p{1});
      end
 
    end % iFunc_McCode constructor
    
    function f = iFunc(self)
      % convert a single iFunc_McCode back to iFunc
      f = iFunc;
      for p = fieldnames(self)'
        f.(p{1}) = self.(p{1});
      end
    end
    
    % overloaded feval which prefers to use 'nan' when axes are undefined
    function [signal, model, ax, name] = feval(self, varargin)
      if numel(varargin) == 0
        [signal, model, ax, name] = feval@iFunc(self, [], nan);
      elseif numel(varargin) == 1 && isvector(varargin{1})
        [signal, model, ax, name] = feval@iFunc(self, varargin{1}, nan);
      else  
        [signal, model, ax, name] = feval@iFunc(self, varargin{:});
      end
    end % feval (override iFunc)
    
    
    % overloaded inputdlg to display the instrument parameters in a dialogue
    % calls private/mccode_run
    function [v,self] = inputdlg(self)
      [v,self] = dialog(self);
    end
    function [v,self] = uitable(self)
      [v,self] = dialog(self);
    end
    
  end
  
end
