function a = permute(a, order)
% PERMUTE Permute object dimensions.
%    B = PERMUTE(A,ORDER) rearranges the dimensions of A so that they
%    are in the order specified by the vector ORDER.  The object produced
%    has the same values as A but the order of the subscripts needed to 
%    access any particular Signal element are rearranged as specified by ORDER.
%    For an N-D object A, numel(ORDER)>=ndims(A). All the elements of 
%    ORDER must be unique. Default ORDER is [ 2 1 ] i.e. transpose.
%
%    PERMUTE on an array of objects permutes the array, not the objects
%    therein.
%
%    PERMUTE is a generalization of transpose (.') 
%
%    You may as well use SQUEEZE to remove singleton dimensions, and
%    RESHAPE to reshape the elements, which total number is kept. To
%    resize the data set, use RESIZE or REDUCEVOLUME.
%
% Example: c=permute(iData(rand(2,3,4)),[2 3 1]);
% Version: $Date$ $Version$ $Author$
% See also iData, iData/size, iData/reshape, iData/resize,
% iData/reducevolume, iData/squeeze

% handle iData array: use built-in permute
if nargin ==1, order=[]; end
if isempty(order), order=[2 1]; end % default is transpose

if numel(a) > 1
  a = builtin(mfilename, a, order);
  return
end

% check if order has the right dimension, else pad with other dimensions
if length(order) < ndims(a)
  for index=1:ndims(a)
    if isempty(find(order==index)), order=[ order index]; end
  end
end

% use permute on Signal, Error, Monitor
if ~isvector(a)
  a = unary(a, 'permute', order);
end

% then swap axes
if length(a.Axes) == length(order)    % all axes defined
  a.Axes = a.Axes(order);
elseif length(a.Axes) % some axes are not defined
  ax = a.Axes;
  for index=(length(ax)+1):length(order)
    ax{index}='';
  end
  a.Axes = ax(order);
end

