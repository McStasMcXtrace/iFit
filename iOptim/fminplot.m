function stop = fminplot(pars, optimValues, state)
% stop=fminplot(pars, optimValues, state) default plotting function showing the criteria evolution 
% as well as main parameters and status

% fields in optimValues: funcount, fval, iteration, procedure
% state values: 'init','interrupt','iter','done'

  persistent parsHistory;
  persistent fvalHistory;
  persistent iterHistory;
  persistent updatePlot
  stop = false;
  
  if ~isfield(optimValues, 'funcount')
    if isfield(optimValues, 'funcCount')
      optimValues.funcount = optimValues.funcCount;
    elseif isfield(optimValues, 'funccount')
      optimValues.funcount = optimValues.funccount;
    end
  end
  
  % check if user has closed the figure to end the optimization
  h = findall(0, 'Tag', 'fminplot');
  if isempty(h) & optimValues.funcount > 5
    stop = true;  % figure was closed: abort optimization by user
    return
  end
  name = [ state ' #' num2str(optimValues.funcount) ' f=' num2str(optimValues.fval,4) ' ' optimValues.procedure ' [close to abort]' ];
  try
    set(h, 'Name', name);
  catch
    stop=true;  % figure is not valid: was closed
    return;
  end
  
  % store data history as rows
  if isempty(fvalHistory) | optimValues.funcount < 1
    parsHistory = pars(:)'; 
    fvalHistory = optimValues.fval;
    iterHistory = 1; % optimValues.funcount;
  else  
    parsHistory = [ parsHistory ; pars(:)' ];
    fvalHistory = [ fvalHistory ; optimValues.fval ];
    iterHistory = [ iterHistory ; length(iterHistory)+1 ]; % [ iterHistory ; optimValues.funcount ]; 
  end
  
  if length(fvalHistory) > 10
    if ~isempty(updatePlot)
      if etime(clock, updatePlot) < 2, return; end % plot every 2 secs
    end
  end
    
  % handle figure
  % only retain one instance of fminplot
  if length(h) > 1, delete(h(2:end)); h=h(1); end
  if isempty(h) & optimValues.funcount <=2 % create it
    h = figure('Tag','fminplot', 'Unit','pixels');
    ishidden = 0;
    tmp = get(h, 'Position'); tmp(3:4) = [500 400];
    set(h, 'Position', tmp);
  end
  
  try
    % raise existing figure (or keep it hidden)
    if gcf ~= h, figure(h); end
  catch
    stop=true;  % figure is not valid: was closed
    return;
  end
  
  % determine best guess up to now
  [dummy, best] = sort(fvalHistory); % sort in ascending order
  best= best(1);
  
  % handle first subplot: criteria
  subplot(1,2,1); % this subplot shows the criteria
  hold off
  g=plot(iterHistory, fvalHistory,'b-', ...
    iterHistory(1),   fvalHistory(1),'ro', ...
    iterHistory(best),fvalHistory(best),'gv', ...
    iterHistory(end), fvalHistory(end), 'rs');
  set(g(end),'MarkerFaceColor','r');
  set(g(end-1),'MarkerFaceColor','g');
  if all(fvalHistory > 0) set(gca, 'yscale', 'log'); end
  xlabel([ 'Nb of Function Evaluations. ' num2str(length(pars)) ' Pars' ]); 
  ylabel('Criteria - {\bf Close figure to abort}');
  if strcmp(state, 'done'),     title('Done','FontWeight','bold','FontSize',14); 
  elseif strcmp(state, 'init'), title('Init','FontWeight','bold','FontSize',14); 
  else                          title('Close figure to abort','FontWeight','bold');  end
  l = legend(g,{'Criteria','Start','Best','Last'},'Location','NorthEast');
  axis auto
  
  % handle second subplot: parameters. The last set is highlighted
  subplot(1,2,2); % this subplot shows some parameters
  hold off
  switch length(pars)
  case 1
    g=plot(fvalHistory,    parsHistory,'bo', ...
        fvalHistory(1),    parsHistory(1),'ro', ...
        fvalHistory(best), parsHistory(best),'gv', ...
        fvalHistory(end),  parsHistory(end),'rs');
    xlabel('FunVal'); ylabel('Par1'); 
  case 2
    g=plot(parsHistory(:,1), parsHistory(:,2),'bo', ...
        parsHistory(1,1),    parsHistory(1,2),'ro', ...
        parsHistory(best,1), parsHistory(best,2),'gv', ...
        parsHistory(end,1),  parsHistory(end,2),'rs');
    xlabel('Par1'); ylabel('Par2'); 
  otherwise
    % find the 3 first parameters that vary
    index = find(sum(abs(diff(parsHistory,1)),1));
    if length(index) > 3, index = index(1:3);
    else
      if ~any(index == 1), index = [ 1 index ]; end
      if ~any(index == 2), index = [ index 2 ]; end
      if ~any(index == 3), index = [ index 3 ]; end
    end
    g=plot3(parsHistory(:,index(1)), parsHistory(:,index(2)),   parsHistory(:,index(3)), 'bo', ...
          parsHistory(1,index(1)),   parsHistory(1,index(2)),   parsHistory(1,index(3)), 'ro', ...
          parsHistory(best,index(1)),parsHistory(best,index(2)),parsHistory(best,index(3)), 'gv', ...
          parsHistory(end,index(1)), parsHistory(end,index(2)), parsHistory(end,index(3)), 'rs');
    xlabel([ 'Par' num2str(index(1)) ]); ylabel([ 'Par' num2str(index(2)) ]); zlabel([ 'Par' num2str(index(3)) ]); 
  end
  
  set(g(end),'MarkerFaceColor','r');
  set(g(end-1),'MarkerFaceColor','g');
  % store userData in case we need to access the optimization history from the plot
  set(gcf,'UserData', struct('parsHistory',parsHistory,'fvalHistory',fvalHistory,'best',best));
  % resize axis frame to display title properly
  p=get(gca,'Position'); p(4)=0.75; set(gca,'Position',p);
  title([ '#' num2str(optimValues.funcount) ' f=' num2str(optimValues.fval,4) ' [close to abort]' sprintf('\n') optimValues.procedure  ]);
  axis auto
  
  updatePlot=clock;
  drawnow

end

