%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <http://en.wikipedia.org/wiki/Absolute_deviation>
function c=least_max(Signal, Error, Model)
% weighted median absolute criteria
% the return value is a scalar
% median(|Signal-Model|/Error)
  c = max(least_absolute(Signal, Error, Model));
end % least_max
