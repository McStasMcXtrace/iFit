function inline_disp(options, funccount, pars, fval)
% function called during minimization procedure
%
% Displays iteration information every 5 steps, then 10 steps, then 100 steps
% or at every step if iteration is negative


  if funccount > 5
    if funccount > 50  & mod(funccount,100) return;
    elseif mod(funccount,10) return; end
  end

  if isfield(options,'Display')
    if strncmp(options.Display, 'iter',4)
      spars=pars(1:min(20,length(pars)));
      spars=mat2str(spars', 4);  % as a row
      if length(spars) > 45, spars=[ spars(1:42) ' ...]' ]; end
      index = isfinite(fval);
      if isfinite(funccount) && ~isempty(index)
        disp(sprintf(' %5.0f    %12.6g          %s', abs(funccount), sum(fval(index)), spars));
      else
        disp(sprintf('                                %s', spars));
      end
    end
  end
  
end % inline_disp
