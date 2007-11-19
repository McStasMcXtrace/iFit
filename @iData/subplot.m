function h=subplot(a, varargin)
% subplot(s) : plot iData array as subplots
%
%   @iData/subplot plot each iData element in a subplot
%
% input:  s: object or array (iData)
% output: h: plot handles (double)
% ex:     subplot([ a a ])
%
% See also iData, iData/plot

% EF 23/11/07 iData implementation

if length(a(:)) == 1
  h=plot(a, varargin);
  return
end

a = squeeze(a); % remove singleton dimensions

if length(size(a)) == 2 & all(size(a) > 1)
  n = size(a,1); m = size(a,2);
else
  p = length(a(:));
  n = floor(sqrt(p));
  m = ceil(p/n);
end

h=[];
for index=1:length(a(:))
  subplot(m,n,index);
  h = [ h plot(a(index)) ];
end
