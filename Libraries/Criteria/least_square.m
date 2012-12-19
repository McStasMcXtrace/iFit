%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <http://en.wikipedia.org/wiki/Least_squares>
function c=least_square(Signal, Error, Model)
% weighted least square criteria, which is also the Chi square
% the return value is a vector, and most optimizers use its sum (except LM).
% (|Signal-Model|/Error).^2
  c = least_absolute(Signal, Error, Model);
  c = c.*c;
end % least_square
