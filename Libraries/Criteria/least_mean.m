function c=least_mean(Signal, Error, Model)
% c=least_mean(Signal, Error, Model)
%
% weighted mean absolute criteria
% the return value is a scalar
% mean(|Signal-Model|/Error)
%
% A good fit corresponds with a criteria lower or equal to 1.
%
% <https://en.wikipedia.org/wiki/Mean_absolute_difference>
% (c) E.Farhi, ILL. License: EUPL.

  c = mean(least_absolute(Signal, Error, Model));
end % least_mean
