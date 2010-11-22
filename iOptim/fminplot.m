function stop = fminplot(pars, optimValues, state)
% default plotting function showing the criteria evolution 
% as well as main parameters and status

% fields in optimValues: funcount, fval, iteration, procedure
% state values: 'init','interrupt','iter','done'

  persistent parsHistory;
  persistent fvalHistory;
  persistent iterHistory;
  stop = false;
  
  % store data history as rows
  if isempty(fvalHistory) | optimValues.iteration <= 1
    parsHistory = pars(:)'; 
    fvalHistory = optimValues.fval;
    iterHistory = optimValues.iteration;
  else  
    parsHistory = [ parsHistory ; pars(:)' ];
    fvalHistory = [ fvalHistory ; optimValues.fval ];
    iterHistory = [ iterHistory ; optimValues.iteration ]; 
  end
  
  if length(fvalHistory) > 9
    if length(fvalHistory) > 999 & mod(length(fvalHistory),1000) return;
    elseif length(fvalHistory) > 99 & mod(length(fvalHistory),100) return;
    elseif mod(length(fvalHistory),5) return; end
  end

    
  % handle figure
  h = findall(0, 'Tag', 'fminplot');
  if length(h) > 1, delete(h(2:end)); h=h(1); end
  if isempty(h) & optimValues.iteration <=2
    h = figure('Tag','fminplot', 'Unit','pixels');
    tmp = get(h, 'Position'); tmp(3:4) = [500 400];
    set(h, 'Position', tmp);
  elseif isempty(h)
    stop = true;  % figure was closed: abort optimization by user
    return
  end
  
  figure(h);
  name = [ localChar(optimValues.procedure) ' ' state ' #' num2str(optimValues.iteration) ];
  try
    set(h, 'Name', name);
  catch
    stop=true;  % figure is not valid: was closed
    return;
  end
  
  % handle first subplot: criteria
  subplot(1,2,1); % this subplot shows the criteria
  g=plot(iterHistory, fvalHistory,'b-', ...
    iterHistory(1), fvalHistory(1),'ro', ...
    iterHistory(end), fvalHistory(end), 'rs');
  set(g(end),'MarkerFaceColor','r');
  xlabel('iteration'); ylabel('criteria'); axis auto
  if strcmp(state, 'done'),     title('Done'); 
  elseif strcmp(state, 'init'), title('Init'); 
  else                          title('Close figure to abort');  end
  
   % handle second subplot: parameters. The last set is highlighted
  subplot(1,2,2); % this subplot shows some parameters
  switch length(pars)
  case 1
    g=plot(fvalHistory, parsHistory,'bo',fvalHistory(end), parsHistory(end),'rs');
    xlabel('FunVal'); ylabel('Par1'); 
  case 2
    g=plot(parsHistory(:,1), parsHistory(:,2),'bo',parsHistory(end,1), parsHistory(end,2),'rs');
    xlabel('Par1'); ylabel('Par2'); 
  otherwise
    g=plot3(parsHistory(:,1), parsHistory(:,2), parsHistory(:,3), 'bo', ...
          parsHistory(end,1), parsHistory(end,2), parsHistory(end,3), 'rs');
    xlabel('Par1'); ylabel('Par2'); ylabel('Par3'); 
  end
  
  set(g(end),'MarkerFaceColor','r');
  title(name);
  axis auto
end
