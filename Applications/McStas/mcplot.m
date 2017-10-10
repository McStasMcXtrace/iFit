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
% (c) E.Farhi, ILL. License: EUPL.

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
  h = subplot(a, 'view2 tight');
  % we identify if this is a scan. Then we add a context menu item to plot 
  % separately each detector
  if all(~cellfun(@isempty, findfield(a, 'xvars'))) && numel(a) == numel(h)
    for index=1:numel(a)
      uicm = get(h(index), 'UIContextMenu');
      if numel(a) == 1, this = a; else this = a(index); end
      [dirname, filename] = fileparts(this.Source);
      component = this.Component;
      if strcmp(component((end-1):end), '_I'), component = component(1:(end-3)); end

      uimenu(uicm, 'Label', [ 'Plot ' component ' scan steps (subplot)...' ], 'Callback', ...
        [ 'tmp_a=iLoad(''' dirname '#' component ''', ''mccode''); ' ...
        'if numel(tmp_a), tmp_a = opensim(iData(tmp_a)); figure; subplot(tmp_a, ''view2 tight''); end; clear tmp_a' ]);
      uimenu(uicm, 'Label', [ 'Stack ' component ' scan steps...' ], 'Callback', ...
        [ 'tmp_a=iLoad(''' dirname '#' component ''', ''mccode''); ' ...
        'if numel(tmp_a), tmp_a = opensim(iData(tmp_a)); tmp_a = cat(0, tmp_a); if ndims(tmp_a) == 3, slice(tmp_a); else figure; plot(tmp_a, ''tight''); end; end; clear tmp_a' ]);
    end
  end

  if ~isempty(options)
    for index=1:numel(options)
      this = options{index};
      saveas(a, '', this);
    end
  end
