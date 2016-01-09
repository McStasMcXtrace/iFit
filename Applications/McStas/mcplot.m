function a=mcplot(varargin)
% mcplot: plot a McCode simulation result (using iFit)
% data = mcplot(filename, ...)
%
% This function displays a McCode simulation result in a single window with 
% subplots. It also returns the resulting objects. 
%
% input:
%  filename: one or more simulation name(s) or directory
%          or a single detector file name
%          if filename does not exist, a file selector is called.
%  optional string: -png, -eps, -fig, -pdf, -jpg 
%            will export further figures directly to files
% 
% output:
%  data: an array of monitors
%
% examples:
%   mcplot
%   mcplot mcstas.sim
%   mcplot -png simulation
%
% SEE ALSO: mcstas, mcplot, mcrun, mcdisplay
% DOC:      Please visit http://www.mcstas.org/

  a = [];
  if nargin ==0
    return
  end
  options={}; % will store output options, such as '-png', '-ps'...
  
  % check if there are output options
  for index=1:numel(varargin)
    this = varargin{index};
    if ~ischar(this), continue; end
    if this(1) == '-'
      options{end+1} = strrep(this, '-', '');
    elseif isdir(this) && ~isempty(dir(fullfile(this,'mccode.sim')))
      varargin{index} = fullfile(this,'mccode.sim');
    end
  end

  % import and display
  a = iData(varargin{:});
  figure;
  subplot(a, 'view2 tight');

  if ~isempty(options)
    for index=1:numel(options)
      this = options{index};
      saveas(a, '', this);
    end
  end
