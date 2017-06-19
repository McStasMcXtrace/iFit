function istop = exec_Callback(pid, Callback, action)
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
      if ~isnan(nout) && nout > 0
        istop = feval(Callback, vars{1:nin});
      elseif ~isnan(nout)
        feval(Callback, vars{1:nin});
      else
        eval(Callback);
      end
    catch ME
      disp(getReport(ME))
    end
  end % Callback
end
