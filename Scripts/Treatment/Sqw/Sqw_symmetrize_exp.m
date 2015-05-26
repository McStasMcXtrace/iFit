function s=Sqw_symmetrize_exp(s)
% function to symmetrise the S(q,w): assume input Sqw contains Bose factor (experimental)
% the return Sqw is classical (symmetric, without T/Bose factor)

  % transform into a symmetric data set in energy axis
  s     = Sqw_deBosify(s);

  % create a new object with an opposite energy axis
  s     = Sqw_symmetrize(s);

