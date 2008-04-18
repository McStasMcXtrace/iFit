function h=subplot(a, varargin)
% subplot(s) : plot iData array as subplots
%
%   @iData/subplot plot each iData element in a subplot
%     subplot(a, [])    uses the best subplot fit
%     subplot(a, [m n]) uses an m x n subplot grid
%
% input:  s: object or array (iData)
%         [m n]: optional subplot grid dimensions
%         additional arguments are passed to the plot method (e.g. color, plot type, ...)
% output: h: plot handles (double)
% ex:     subplot([ a a ])
%
% See also iData, iData/plot

% EF 23/11/07 iData implementation

if length(a(:)) == 1
  h=plot(a, varargin{:});
  return
end

a = squeeze(a); % remove singleton dimensions
m=[];
n=[];
if length(varargin) >=1
  if isnumeric(varargin{1}) | isempty(varargin{1})
    dim = varargin{1};
    if isempty(dim)
      % will use best fit
    elseif length(dim) == 1 & dim(1) > 0
      m = dim; 
    else m=dim(1); n=dim(2); end
    if length(varargin) >= 2  
      varargin = varargin(2:end);
    else varargin = {}; end
  elseif length(size(a)) == 2 & any(size(a) > 1)
    m = size(a,1); n = size(a,2);
  end
end
if any(m==0), m=[]; end
if isempty(m)
  p = length(a(:));
  n = floor(sqrt(p));
  m = ceil(p/n);
elseif isempty(n)
  n = ceil(length(a(:))/m);
end

h=[];
for index=1:length(a(:))
  if ~isempty(a(index))
    subplot(m,n,index);
    h = [ h plot(a(index), varargin{:}) ];
  else h = [ h nan ];
  end
end
