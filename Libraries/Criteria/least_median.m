%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <http://en.wikipedia.org/wiki/Median_absolute_deviation>
function c=least_median(Signal, Error, Model)
% weighted median absolute criteria
% the return value is a scalar
% median(|Signal-Model|/Error)
  c = median(least_absolute(Signal, Error, Model));
end % least_median
