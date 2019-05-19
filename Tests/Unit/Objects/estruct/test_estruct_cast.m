function ret=test_estruct_cast
  s = estruct(rand(5));
  c = cast(s, 'double');
  if isnumeric(c), ret = [ 'OK ' mfilename ];
  else ret = 'FAILED'; end
