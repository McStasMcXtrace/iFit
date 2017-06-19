function istop = exec_Callback(pid, Callback, action)
  % ececutes a Callback in a reduced environment.
  istop = 0; % failed ExternalFcns ignored
  if ~isempty(Callback) && (ischar(Callback) || isa(Callback, 'function_handle'))
    if isa(Callback, 'function_handle')
      nin = nargin(Callback);
      nout= nargout(Callback);
    else
      nin = 0; nout = nan;
    end
    
    if nin
      vars = { pid, action };
    end
    try
      if ~isnan(nout) && nout
        try
          istop = feval(Callback, vars{1:nin});
        catch
          feval(Callback, vars{1:nin});
        end
      elseif ~isnan(nout)
        feval(Callback, vars{1:nin});
      else
        try
          istop = eval(Callback);
        catch
          eval(Callback);
        end
      end
    catch ME
      disp(getReport(ME))
    end
    if istop
      exit_Process(pid, 'kill');
    end
  end % Callback
end
