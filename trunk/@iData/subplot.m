function h=subplot(a, varargin)
% h = subplot(s) : plot iData array as subplots
%
%   @iData/subplot plot each iData element in a subplot
%     subplot(a, [])    uses the best subplot fit
%     subplot(a, [m n]) uses an m x n subplot grid
%     subplot(a, [m n], options) sends options to the plot
%
% input:  s: object or array (iData)
%         [m n]: optional subplot grid dimensions
%         additional arguments are passed to the plot method (e.g. color, plot type, ...)
% output: h: plot handles (double)
% ex:     subplot([ a a ])
%
% Version: $Revision: 1.15 $
% See also iData, iData/plot

% EF 23/11/07 iData implementation

a = squeeze(a); % remove singleton dimensions
m=[]; n=[]; dim=[];
% handle optional arguments
method = '';
for index=1:length(varargin)
  if ischar(varargin{index}),        method = varargin{index};
  elseif isa(varargin{index},'iData') 
    a = [a(:) ; varargin{index} ];
  elseif isnumeric(varargin{index}), dim    = varargin{index};
  end
end
if length(a(:)) == 1
  h=plot(a, method);
  return
end

if length(dim) == 1 & dim(1) > 0
  m = dim; 
elseif length(dim) == 2, m=dim(1); n=dim(2); 
else m=[]; end
  
if any(m==0), m=[]; end
if isempty(m)
  if length(size(a)) == 2 & all(size(a) > 1)
    m = size(a,1); n = size(a,2);
  else
    p = length(a(:));
    n = floor(sqrt(p));
    m = ceil(p/n);
  end
elseif isempty(n)
  n = ceil(length(a(:))/m);
end

h=[];
for index=1:length(a(:))
  if ~isempty(a(index))
    subplot(m,n,index);
    if length(a(:)) > 12, 
      % compact layout
      method =[ method ' hide_axes' ]; 
    end
     
    this_h = plot(a(index), method);
    
    h = [ h this_h ];
  else h = [ h nan ];
  end
end
