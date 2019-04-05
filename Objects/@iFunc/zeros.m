function s = zeros(iFunc_in,varargin)
% s = zeros(s,N,M,P,...) : initialize an iFunc array
%
%   @iFunc/zeros function to create an array of 's' iFunc objects
%   The object 's' is duplicated into an array. Use s=iFunc to get an empty array.
%
% input:  s: object or array (iFunc)
% output: b: object or array (iFunc)
% ex: zeros(iFunc,5,5) will create a 5-by-5 empty iFunc array
%     zeros(s,5,5) will return a 5-by-5 array filled with 's'
%
% Version: $Date$ $Version$ $Author$
% See also iFunc

% EF 27/07/00 creation
% EF 23/09/07 iFunc implementation

if nargin == 1, 
  if numel(iFunc_in) == 1, s=iFunc; return; end
  s=zeros(size(iFunc_in));
  iFunc_in = iFunc;
else s = zeros(varargin{:}); end

if isempty(s), s=iFunc; return; end

s_index = 1;
long_s  = iFunc_in(1);

for i = 2:numel(s)
  s_index = s_index +1;
  if s_index > numel(iFunc_in), s_index = 1; end
  long_s = [ long_s copyobj(iFunc_in(s_index)) ];
end

s = reshape(long_s, size(s));

