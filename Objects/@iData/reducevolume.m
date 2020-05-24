function s = reducevolume(a, R)
% REDUCEVOLUME reduce an object size
%   B = REDUCEVOLUME(A, [Rx Ry Rz ... ]) reduces the number
%     of elements in the object by only keeping every Rx element in the x
%     direction, every Ry element in the y direction, and every Rz element
%     in the z direction and so on.
%
%   B = REDUCEVOLUME(A, R)
%     If a scalar R is used to indicate the amount or
%     reduction instead of a vector, the reduction is assumed to
%     be R on all axes. 
%
%   B = REDUCEVOLUME(A)
%     When omitted, the volume/size reduction is performed on bigger axes until
%     the final object contains less than 1e6 elements.
%
% You may also use SQUEEZE to remove singleton dimensions, RESIZE to
% compress/expand its size, and PACK to compress its storage.
%
% Example: a=iData(peaks(2000)); b=reducevolume(a); prod(size(b)) < 2e6
% Version: $Date$ $Version$ $Author$
% See also iData, iData/squeeze, iData/pack, iData/resize, iData/size

% rebinning object so that its number of elements is smaller 

if nargin == 1
  R = []; % will guess to reduce down to 1e6
end

% handle input iData arrays
if numel(a) > 1
  s = [];
  for index=1:numel(a)
    s =  [ s feval(mfilename, a(index), R) ];
  end
  s = reshape(s, size(a));
  return
end

% determine best reduction factor
if ischar(R) || isempty(R)
  S  = size(a);
  R  = ones(size(S));
  S0 = S;
  
  % loop until we get les than 1e6 elements
  while prod(S) > 1e6
    % identify the biggest axis, and reduce it by increasing R
    for index=1:length(R)
      [dummy, j] = sort(S);
      % S(j(end)) is the biggest element
      R(j(end)) = R(j(end))+1;
      S = S0 ./ R;
    end
  end
 
end

% case of event data sets
if isvector(a)
  if length(R) > 1, R = max(R); end
elseif isnumeric(R) && length(R) == 1
  % R is scalar and object is not vectorial
  R = R*ones(1, ndims(a));
end

S = [];
S.type='()';
S.subs=cell(1,length(R));

% scan dimensions and rebin them
myisvector = @(c)length(c) == numel(c);
for index=1:length(R)
  x = getaxis(a, index);
  if myisvector(x), lx=length(x);
  else              lx=size(x,index); end
  if R(index) > 1
    S.subs{index} = ceil(1:R(index):lx);
  else
    S.subs{index} = ':';
  end
end
clear x

% get sub-object
s = subsref(a, S);
