function y=size(s, varargin)
% size(s) : get estruct object size (number of elements)
%
%   @estruct/size function to get estruct object size
%
% input:  s: object or array (estruct)
%         dim: optional dimension/rank to inquire
% output: v: size of the estruct Signal (double)
% ex:     size(estruct), size(estruct,1)
%
% Version: $Date$
% See also estruct, estruct/get, length

if numel(s) > 1  % this is an array of estruct
  y = builtin('size', s, varargin{:});
  return
end

if numel(s) == 0
  y=[0 0];
else
  % use cache when available for faster execution
  if isfield(s.Private,'cache') && isfield(s.Private.cache,'size') && ~isempty(s.Private.cache.size)
    y = s.Private.cache.size;
  else
    y = size(subsref(s,struct('type','.','subs','Signal')));
  end
  if nargin > 1 && all(cellfun(@isnumeric, varargin))
    y = y(varargin{:}); 
  end
end

