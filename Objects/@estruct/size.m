function y=size(s, varargin)
% SIZE Get object size.
%   D = SIZE(X) for a single object returns the size of its Signal. For
%   an array of objects, the size of the array is returned.
%
%   M = SIZE(X,DIM) returns the length of the dimension specified
%   by the scalar DIM.  For example, SIZE(X,1) returns the number
%   of rows. If DIM > NDIMS(X), M will be 1.
%
%   To get the size of all objects in an array S, use ARRAYFUN('size',S)
%
% Example: s=estruct(rand(5)); all(size(s) == size(s.Signal))
% Version: $Date$ $Version$ $Author$
% See also estruct, get, length, ndims

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
    s.Private.cache.size = y;
  end
  try
    y = y(varargin{:});
  catch
    if s.verbose
      warning([ mfilename ': ERROR: Invalid dimension specification [ ' sprintf('%i ', varargin{:}) '] for object ' s.Tag ' ' s.Name ' of dimensionality ' num2str(ndims(s)) ])
    end
    y = [];
  end
end
