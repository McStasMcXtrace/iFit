classdef iFunc_McCode < iFunc

  % iFunc_McCode: create an iFunc_McCode from e.g. an instrument model
  %
  % The iFunc_McCode class is a model holding a ray-tracing instrument simulation.
  %
  % Example: s=mccode('defaults')
  %
  % Useful methods for this iFunc flavour:
  %
  % methods(iFunc_McCode)
  %   all iFunc methods can be used.
  % iFunc_McCode(a)
  %   convert a 4D model [a=iFunc class] into an iFunc_McCode to give access to
  %   the methods below.
  % plot(a)
  %   plot the instrument geometry
  % subplot(a)
  %   plot the instrument simulation results (monitors)
  % publish(a)
  %   Generate a readable document with all results
  % edit(a)
  %   Show the instrument code, and re-compile the instrument when modified
  % dialog(a)
  %   Display a simple dialogue window to run a single or series of simulations
  %
  % input:
  %   can be an iFunc or struct or any set of parameters to generate a McCode object.
  %   when not given an iFunc, the parameters to mccode are expected.
  %
  % output: an iFunc_McCode object
  %
  % See also: mccode, iFunc, iFunc_McCode/plot, iFunc_McCode/subplot, iFunc_McCode/dialog, iFunc_McCode/edit, iFunc_McCode/publish

  properties
  end

  methods
    function obj = iFunc_McCode(varargin)
      % iFunc_McCode: create the iFunc_McCode subclass
      %
      % input:
      %   can be an iFunc or any set of parameters to generate a McCode object
      obj = obj@iFunc;
      if nargin == 0
        m = mccode('defaults');
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
        obj = m; obj.class = mfilename;
        return;
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
      obj.class = mfilename;
 
    end % iFunc_McCode constructor
    
    function f = iFunc(self)
      % iFunc_McCode: iFunc: convert a single iFunc_McCode back to iFunc
      f = iFunc;
      for p = fieldnames(self)'
        f.(p{1}) = self.(p{1});
      end
      f.class = 'iFunc';
    end
    
    % overloaded feval which prefers to use 'nan' when axes are undefined
    function [signal, self, ax, name] = feval(self, varargin)
      % iFunc_McCode: feval: Evaluate an McCode model value using given parameters and axes.
      % When no axis is given, the raw monitor values are returned.
      % To get the result as iData objects, use: iData(model, ...)
      %   or get model.UserData.monitors after the calculation.
      if numel(varargin) == 0
        [signal, self, ax, name] = feval@iFunc(self, [], nan);
      elseif numel(varargin) == 1 && isvector(varargin{1})
        [signal, self, ax, name] = feval@iFunc(self, varargin{1}, nan);
      else  
        [signal, self, ax, name] = feval@iFunc(self, varargin{:});
      end
      if ~isempty(inputname(1))
        assignin('caller',inputname(1),self); % update in original object
      end
    end % feval (override iFunc)
    
    
    % overloaded inputdlg to display the instrument parameters in a dialogue
    % calls private/mccode_run
    function [v,self] = inputdlg(self)
    % iFunc_McCode: inputdlg: run a McCode instrument model and request instrument parameters in a dialogue
      [v,self] = dialog(self);
    end
    function [v,self] = uitable(self)
    % iFunc_McCode: uitable: run a McCode instrument model and request instrument parameters in a dialogue
      [v,self] = dialog(self);
    end
    
    % overloaded edit to edit the instrument code
    function self = edit(self)
      % iFunc_McCode: edit: edit the instrument source description, and re-build the model when modified.
      
      % we display the instrument source code (when not 'out') in TextEdit
      % as a temporary file. We wait for TextEdit to close.
      
      source = self.UserData.instrument_source;
      if ~ischar(source), return; end
      
      d = tempname;
      mkdir(d);
      filename = fullfile(d, self.UserData.options.instrument);
      % write the comtent of the current object in the temporary file
      fid=fopen(filename, 'w');
      if fid==-1, 
        disp([ mfilename ': WARNING: model ' self.Name ' ' self.Tag ' could not write ' filename ]);
        return
      else 
        fprintf(fid, '%s\n', source);
        fclose(fid);
      end
      
      fig = TextEdit(filename);
      waitfor(fig);
      % get the (updated) source file and create new object
      source2 = fileread(filename);
      % test if the new source has been updated
      n = numel(source);
      has_changed = false;
      if ~all(source(1:n) == source2(1:n))  % common part
        has_changed = true;
      end
      if numel(source2) > n && ~all(isspace(source2((n+1):end)))
        has_changed = true;
      end
      
      if has_changed
        % compile new object
        self = iFunc_McCode(filename);
        if nargin == 0 && ~isempty(inputname(1))
          assignin('caller',inputname(1),self);
        end
      end
    end % edit
    
    % overloaded subplot
    function h = subplot(self, varargin)
      if numel(self) == 1
        h = subplot(self.UserData.monitors, varargin{:});
      else
        h = subplot@iFunc(self, varargin{:});
      end
    end % subplot
    
  end
  
end
