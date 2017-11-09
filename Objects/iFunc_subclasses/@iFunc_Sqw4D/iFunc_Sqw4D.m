classdef iFunc_Sqw4D < iFunc

  properties
  end
  
  methods
  
    function obj = iFunc_Sqw4D(varargin)
      % create an iFunc_Sqw4D from e.g. an iFunc 4D object
      %
      % input:
      %   can be an iFunc or any set of parameters to generate a Sqw4D object.
      %   when not given an iFunc, the parameters to sqw_phonons are expected.
      %
      % See also: sqw_phonons, sqw_cubic_monoatomic
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
    end % iFunc_Sqw4D constructor
    
    function [fig, s, k]=plot(self, varargin)
      % iFunc_Sqw4D: plot: plot dispersions along principal axes and vDOS
      if isempty(varargin) || (numel(varargin) == 1 && ischar(varargin{1}))
        if isempty(varargin), varargin{1} = 'plot meV'; end
        [s,k,fig]=sqw_kpath(self, varargin{1});
      else
        fig = figure;
        fig = plot@iFunc(self, varargin{:});
      end
      if ~isempty(inputname(1))
        assignin('caller',inputname(1),self); % update in original object
      end
    end % plot
    
    function h=plot3(s, varargin)
      % iFunc_Sqw4D: plot3: plot a 3D view of the dispersions in H=0 plane
      s = maxfreq(s);
      % evaluate the 4D model onto a mesh filling the Brillouin zone [-0.5:0.5 ]
      s.UserData.DOS     = [];  % make sure we re-evaluate again on a finer grid
      s.UserData.maxFreq = max(s.UserData.maxFreq(:));
      qk=linspace(0,0.5,50); qh=0; ql=qk; 
      w =linspace(0.01,s.UserData.maxFreq*1.2,51);
      f =iData(s,[],qh,qk,ql',w);
      % plot in 3D
      h = plot3(log(f(1,:, :,:)), varargin{:}); % h=0
      if ~isempty(inputname(1))
        assignin('caller',inputname(1),s); % update in original object
      end
    end % plot3
    
    function h=scatter3(s, varargin)
      % iFunc_Sqw4D: scatter3: plot a 3D scatter view of the dispersions in H=0 plane
      s = maxfreq(s);
      % evaluate the 4D model onto a mesh filling the Brillouin zone [-0.5:0.5 ]
      s.UserData.DOS     = [];  % make sure we re-evaluate again on a finer grid
      s.UserData.maxFreq = max(s.UserData.maxFreq(:));
      qk=linspace(0,0.5,50); qh=0; ql=qk; 
      w =linspace(0.01,s.UserData.maxFreq*1.2,51);
      f =iData(s,[],qh,qk,ql',w);
      % plot in 3D
      h = scatter3(log(f(1,:, :,:)), varargin{:}); % h=0
      if ~isempty(inputname(1))
        assignin('caller',inputname(1),s); % update in original object
      end
    end % scatter3
  
    % methods for Sqw 4D
    
      
    % kpath
    % thermochemistry   
    % bosify
    % debosify
    % gdos
    % powder
    % sq
    % publish -> report
  end % methods
  
end % classdef
  

% private functions used in the class ----------------------------------------

function s = maxfreq(s)
  % get a quick estimate of the max frequency
  if  ~isfield(s.UserData,'maxFreq') || isempty(s.UserData.maxFreq) ...
    || all(s.UserData.maxFreq <= 0)
    qh=linspace(-.5,.5,10);qk=qh; ql=qh; w=linspace(0.01,50,11);
    f=iData(s,[],qh,qk,ql',w);
    if isfield(s.UserData, 'FREQ') && ~isempty(s.UserData.FREQ)
      s.UserData.maxFreq = max(s.UserData.FREQ(:));
      disp([ mfilename ': maximum phonon energy ' num2str(max(s.UserData.maxFreq)) ' [meV] in ' s.Name ]);
    end
    if ~isfield(s.UserData, 'maxFreq') || isempty(s.UserData.maxFreq) ...
      || ~isfinite(s.UserData.maxFreq) || s.UserData.maxFreq <= 0
      s.UserData.maxFreq = 100;
    end
    
  end
end

