function c = mtimes(a,b)
% c = mtimes(a,b) : computes the matrix product of estruct objects
%
%   @estruct/mtimes function to compute the matrix product of data sets=a*b
%      the number of columns of a must equal the number of rows of b.

%
% input:  a: object or array (estruct or numeric)
%         b: object or array (estruct or numeric)
% output: c: object or array (estruct)
% ex:     c=a*2;
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/minus, estruct/plus, estruct/times, estruct/rdivide, estruct/power

if nargin == 1,
  b = a;
end
c=[];

% handle handle array as input
if numel(b) > 1
  c = zeros(estruct, length(b), 1);
  for index=1:length(b)
    c(index) = feval(mfilename, a, b(index));
  end
  return
elseif numel(a) > 1
  c = zeros(estruct, length(a), 1);
  for index=1:length(a)
    c(index) = feval(mfilename, a(index), b);
  end
  return
end

if isscalar(a) | isscalar(b)
  c = binary(a, b, 'times');
elseif ndims(a) == 2 & ndims(b) == 2
  if isa(a, 'estruct') && numel(a) > 1, a=a(1); end
  if isa(b, 'estruct') && numel(b) > 1, b=b(end); end

  if size(a,2) ~= size(b,1)
    error([ mfilename ': the number of columns of a (' num2str(size(a,2)) ') must equal the number of rows of b (' num2str(size(b,1)) ').' ]);
  end
  cmd=a.Command;
  c = copyobj(a);
  c1 = dot(getaxis(a,2) , getaxis(b,1));
  c2 = dot(getaxis(a,1) , getaxis(b,2));
  c0 = double(a)  * double(b);
  [a1n] = getaxis(a, '1'); a1l = label(a,1);
  [a2n] = getaxis(a, '2'); a2l = label(a,2);
  [a0n] = getaxis(a, '0'); a0l = label(a,0);
  [b1n] = getaxis(b, '1'); b1l = label(b,1);
  [b2n] = getaxis(b, '2'); b2l = label(b,2);
  [b0n] = getaxis(b, '0'); b0l = label(b,0);
  ae = get(a, 'Error'); am = get(a, 'Monitor');
  be = get(b, 'Error'); bm = get(b, 'Monitor');
  setalias(c,[ a2n '_' b1n ],c1); label(c,[ a2n '_' b1n ],[ a2l ' ' b1l ]);
  setalias(c,[ a1n '_' b2n ],c2); label(c,[ a1n '_' b2n ],[ a1l ' ' b2l ]);
  setalias(c, 'Signal', c0); label(c, 'Signal', [ a0l '*' b0l ]);
  setalias(c, 'Error',   ae * be);
  setalias(c, 'Monitor', am * bm);
  c.Command=cmd;
  c = history(c, mfilename, a,b);
else
  error([ mfilename ': Matrix estruct multiplication only supported for matrix estruct object (ndims=2), and not current ndims(a)=' num2str(ndims(a)) ' and ndims(b)=' num2str(ndims(b)) '.' ]);
end


