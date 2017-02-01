function stop = fminstop(x, optimValues, state)
% fminstop: an OutputFcn / PlotFcns function which displays a 'STOP' button
%
% this later can be pressed during a fit to abort the current fit procedure.
%
% Other defined function to be used as OutputFcn and PlotFcns:
%   @fminplot           shows criteria and parameter space, with a STOP button
%   @optimplotx         plots the current point
%   @optimplotfval      plots the function value
%   @optimplotfunccount plots the function count
%
% example: 
%   fmin(@objective, [], 'OutputFcn=fminplot')

persistent fig

stop = false;

switch state
case 'init'
  % open a figure with a single toggle button
  fig = findall(0, 'Tag','Optim:fminstop');
  if length(fig) > 1, delete(fig(2:end)); end % unique instance

  if isempty(fig) % create window in case state was never 'init'
    fig = fminstop_create;
  end
  
case {'iter','interrupt'}
  % return stop if the figure is closed
  if ~isempty(fig) && ~ishandle(fig) % closed by user
    stop = true;
    fig  = [];
    return; 
  end

  if isempty(fig) % create window in case state was never 'init'
    fig = fminstop_create;
  end
case 'done'
  % close window
  close(fig);
  stop = true;
  fig = [];
end

% ------------------------------------------------------------------------------
function fig = fminstop_create
  fig = figure('Tag','Optim:fminstop','MenuBar','None','NextPlot','new', ...
      'HandleVisibility','callback','Name','Fit [close to abort]');
  p = get(fig, 'Position');
  p(3:4) = [100 50];
  set(fig, 'Position',p);
  h = uicontrol(fig,'String','STOP FIT',...
    'Style','pushbutton','callback','fminstop([],[],''done'');','BackgroundColor','red', ...
    'Units','normalized','Position',[.1 .1 .8 .8], ...
    'ToolTipString','Click me to abort current fit');
    
  drawnow;
