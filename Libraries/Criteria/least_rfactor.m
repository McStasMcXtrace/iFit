function c=least_rfactor(Signal, Error, Model)
% c=least_rfactor(Signal, Error, Model)
%
% minimal R-factor (weighted least square criteria) normalised to the Signal.
% the return value is a vector, and most optimizers use its sum (except LM).
% (|Signal-Model|/Error).^2/(Signal/Error).^2
%
% A good fit corresponds with an R-factor lower than .2
%
% <http://en.wikipedia.org/wiki/R-factor_%28crystallography%29>
% (c) E.Farhi, ILL. License: EUPL.

  if ~isnumeric(Signal) || ~isnumeric(Model), return; end
  if isempty(Error) || isscalar(Error) || all(Error == Error(end))
    index = find(isfinite(Model) & isfinite(Signal) & Signal);
    c = sqrt(((Signal(index)-Model(index)).^2)./(Signal(index).^2)); % raw normalised least absolute
  else
    index = find(isfinite(Error) & isfinite(Model) & isfinite(Signal) & Signal);
    residuals  = Signal - Model;
    % make sure weight=1/sigma does not reach unrealistic values
    %   initially, most weights will be equal, but when fit impproves, 
    %   stdE will get lower, allowing better matching of initial weight.
    normE = sum(Error(index));
    stdE  = std(residuals(index));
    Error( Error < stdE/4 ) = stdE/4; 
    Error = Error *(normE/sum(Error(index)));
    
    if isempty(index), c=Inf;
    else               c=sqrt(((residuals(index)./Error(index)).^2)./((Signal(index)./Error(index)).^2));
    end
  end
end % least_rfactor
