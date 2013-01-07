function c=least_absolute(Signal, Error, Model)
% weighted least absolute criteria
% the return value is a vector, and most optimizers use its sum (except LM).
% |Signal-Model|/Error
%
% A good fit corresponds with a criteria lower or equal to 1.
%
% <http://en.wikipedia.org/wiki/Least_absolute_deviation>

  if isempty(Error) || isscalar(Error) || all(Error == Error(end))
    index = find(isfinite(Model) & isfinite(Signal));
    c = abs(Signal(index)-Model(index)); % raw least absolute
  else
    % find minimal non zero Error
    Error = abs(Error);
    index = find(Error~=0 & isfinite(Error));
    minError = min(Error(index));
    % find zero Error, which should be replaced by minimal Error whenever possible
    index = find(Error == 0);
    Error(index) = minError;
    index = find(isfinite(Error) & isfinite(Model) & isfinite(Signal));
    if isempty(index), c=Inf;
    else               c=abs((Signal(index)-Model(index))./Error(index));
    end
  end
end % least_absolute
