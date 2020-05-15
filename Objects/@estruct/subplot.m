function h=subplot(a, varargin)
% SUBPLOT  Display objects as a grid in tiled positions.
%   SUBPLOT(A) and SUBPLOT(a, []) uses the best subplot fit. The object A can be
%   given as an array, such as in SUBPLOT([A B C...]), or as separate argumentys
%   SUBPLOT(A,B,C...).
%
%   SUBPLOT(A, [m n]) uses an m x n subplot grid.
%
%   SUBPLOT(A, [m n], options) sends options to the plot.
%
% Example: a=estruct(peaks); h=subplot([a a]); tf=ishandle(h); delete(h); all(tf)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plot

a = squeeze(a); % remove singleton dimensions
m=[]; n=[]; dim=[];
% handle optional arguments
method = '';
varg = {};
for index=1:length(varargin)
  if ischar(varargin{index}),        method = varargin{index};
  elseif isa(varargin{index},class(a)) 
    if numel(a) == 1
      a = [a ; varargin{index} ];
    else
      a = [a(:) ; varargin{index} ];
    end
    varargin{index} = [];
  elseif isempty(dim) && isnumeric(varargin{index}) && numel(varargin{index}) == 2
    dim    = varargin{index};
  else varg{end+1} = varargin{index};
  end
end
clear varargin
if numel(a) == 1
  h=plot(a, method, varg{:});
  return
end

if length(dim) == 1 && dim(1) > 0
  m = dim;
elseif length(dim) == 2, m=dim(1); n=dim(2); 
else m=[]; end
  
if any(m==0), m=[]; end
if isempty(m)
  if length(size(a)) == 2 & all(size(a) > 1)
    m = size(a,1); n = size(a,2);
  else
    p = numel(a);
    n = floor(sqrt(p));
    m = ceil(p/n);
  end
elseif isempty(n)
  n = ceil(length(a(:))/m);
end

h=[];
for index=1:numel(a)
  if ~isempty(a(index))
    if numel(a) > 12, 
      % compact layout
      method =[ method ' hide_axes' ]; 
    end
    subplot(m,n,index);
     
    this_h = plot(a(index), method, varg{:});
    
    h = [ h ; this_h(:) ];
  else h = [ h ; nan ];
  end
end
