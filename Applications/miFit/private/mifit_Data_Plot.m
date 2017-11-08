function mifit_Data_Plot(style, varargin)
% Data/Sub Plot (Overview)
  % get the Data from either 'style' as 1st arg or further varargin
  if    ~nargin, style = ''; end
  if     nargin     && isa(style, 'iData'),       d=style;
    if nargin > 1   && ischar(varargin{1}),   style=varargin{1}; 
    else style = ''; end
  elseif nargin > 1 && isa(varargin{1}, 'iData'), d=varargin{1};
  else d = mifit_List_Data_pull; end
  if iscell(d), d= [ d{:} ]; end
  if ~nargin || isempty(style) || ~ischar(style), style=''; end
  if all(isempty(d)), return; end
  set(mifit_fig,'Pointer','watch');
  f=figure('Visible','off');
  subplot(d,[ style ' grid tight replace' ]);
  if f~=gcf, close(f); else set(f,'Visible','on'); end
  set(mifit_fig,'Pointer','arrow');
