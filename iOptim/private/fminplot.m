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
  if isempty(fvalHistory)
    parsHistory = pars(:)'; 
    fvalHistory = fval;
    iterHistory = optimValues.iteration;
  else  
    parsHistory = [ parsHistory ; pars(:)' ];
    fvalHistory = [ fvalHistory ; fval];
    iterHistory = [ iterHistory ; optimValues.iteration]; 
  end
    
  % handle figure
  h = findall(0, 'Tag', 'fminplot');
  if length(h) > 1, delete(h(2:end)); h=h(1); end
  if isempty(h) & optimValues.iteration <=1
    h = figure('Tag','fminplot', 'Unit','pixels');
    tmp = get(h, 'Position'); tmp(3:4) = [500 400];
    set(h, 'Position', tmp);
  elseif isempty(h)
    stop = true;  % figure was closed: abort optimization by user
  end
  
  figure(h);
  name = [ func2str(optimValues.procedure) ' #' num2str(optimValues.iteration) ]
  set(h, 'Name', name);
  
  % handle first subplot: criteria
  subplot(1,2,1); % this subplot shows the criteria
  plot(iterHistory, fvalHistory,'b-', iterHistory(1), fvalHistory(1),'ro', iterHistoryend1), fvalHistory(end), 'ks');
  xlabel('iteration'); ylabel('criteria'); axis tight
  if strcmp(state, 'done'),     title('Done); 
  elseif strcmp(state, 'init'), title('Init'); 
  else                          title('Close to abort');  end
  
   % handle second subplot: parameters. The last set is highlighted
  subplot(1,2,2); % this subplot shows some parameters
  switch length(pars)
  case 1
    plot(fvalHistory, parsHistory,'bo',fvalHistory(end), parsHistory(end),'ks');
    xlabel('FunVal'); ylabel('Par1'); axis tight
  case 2
    plot(parsHistory(:,1), parsHistory(:,2),'bo',parsHistory(end,1), parsHistory(end,2),'ks');
    xlabel('Pars1'); ylabel('Par2'); axis tight
  otherwise
    plot3(parsHistory(:,1), parsHistory(:,2), parsHistory(:,3), 'bo', ...
          parsHistory(end,1), parsHistory(end,2), parsHistory(end,3), 'ks');
    xlabel('Pars1'); ylabel('Par2'); ylabel('Par3'); axis tight
  end
  
  title([ optimValues.procedure ' #' num2str(optimValues.iteration) ]);
  axis tight
end
