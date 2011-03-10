function stop = fminplot(pars, optimValues, state)
% default plotting function showing the criteria evolution 
% as well as main parameters and status

% fields in optimValues: funcount, fval, iteration, procedure
% state values: 'init','interrupt','iter','done'

  persistent parsHistory;
  persistent fvalHistory;
  persistent iterHistory;
  stop = false;
  
  % check if user has closed the figure to end the optimization
  h = findall(0, 'Tag', 'fminplot');
  if isempty(h) & optimValues.funcount > 5
    stop = true;  % figure was closed: abort optimization by user
    return
  end
  
  % store data history as rows
  if isempty(fvalHistory) | optimValues.funcount <= 1
    parsHistory = pars(:)'; 
    fvalHistory = optimValues.fval;
    iterHistory = optimValues.funcount;
  else  
    parsHistory = [ parsHistory ; pars(:)' ];
    fvalHistory = [ fvalHistory ; optimValues.fval ];
    iterHistory = [ iterHistory ; optimValues.funcount ]; 
  end
  
  if length(fvalHistory) > 9
    if length(fvalHistory) > 999 & mod(length(fvalHistory),1000) return;
    elseif length(fvalHistory) > 99 & mod(length(fvalHistory),100) return;
    elseif mod(length(fvalHistory),5) return; end
  end

    
  % handle figure
  if length(h) > 1, delete(h(2:end)); h=h(1); end
  if isempty(h) & optimValues.funcount <=2
    h = figure('Tag','fminplot', 'Unit','pixels');
    tmp = get(h, 'Position'); tmp(3:4) = [500 400];
    set(h, 'Position', tmp);
  end
  
  figure(h);
  name = [ optimValues.procedure ' ' state ' #' num2str(optimValues.funcount) ' [close to abort]' ];
  try
    set(h, 'Name', name);
  catch
    stop=true;  % figure is not valid: was closed
    return;
  end
  
  % determine best guess up to now
  [dummy, best] = sort(fvalHistory); % sort in ascending order
  best= best(1);
  
  % handle first subplot: criteria
  subplot(1,2,1); % this subplot shows the criteria
  g=plot(iterHistory, fvalHistory,'b-', ...
    iterHistory(1), fvalHistory(1),'ro', ...
    iterHistory(best),fvalHistory(best),'gv', ...
    iterHistory(end), fvalHistory(end), 'rs');
  set(g(end),'MarkerFaceColor','r');
  set(g(end-1),'MarkerFaceColor','g');
  if all(fvalHistory > 0) set(gca, 'yscale', 'log'); end
  xlabel('Nb of Function Evaluations'); ylabel('Criteria - {\bf Close} figure to abort'); axis auto
  if strcmp(state, 'done'),     title('Done','FontWeight','bold','FontSize',14); 
  elseif strcmp(state, 'init'), title('Init','FontWeight','bold','FontSize',14); 
  else                          title('Close figure to abort','FontWeight','bold');  end
  l = legend(g,{'Criteria','Start','Best','Last'},'Location','NorthEast');
  
   % handle second subplot: parameters. The last set is highlighted
  subplot(1,2,2); % this subplot shows some parameters
  switch length(pars)
  case 1
    g=plot(fvalHistory,    parsHistory,'bo', ...
        fvalHistory(1), parsHistory(1),'ro', ...
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
    g=plot3(parsHistory(:,1), parsHistory(:,2), parsHistory(:,3), 'bo', ...
          parsHistory(1,1), parsHistory(1,2), parsHistory(1,3), 'ro', ...
          parsHistory(best,1), parsHistory(best,2), parsHistory(best,3), 'gv', ...
          parsHistory(end,1), parsHistory(end,2), parsHistory(end,3), 'rs');
    xlabel('Par1'); ylabel('Par2'); ylabel('Par3'); 
  end
  
  set(g(end),'MarkerFaceColor','r');
  set(g(end-1),'MarkerFaceColor','g');
  title(name);
  axis auto
end
