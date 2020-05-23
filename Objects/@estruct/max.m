function [m,id] = max(a,b, dim)
% MAX    Largest component.
%   MAX(A) returns a single value as the maximum value of the object signal.
%
%   MAX(A,B) returns an object which signal is the highest of A and B.
%
%   MAX(A,[], DIM) returns maximum value along dimension 'DIM'
%
%   [M,I] = MAX(A,...) returns the maximum value/object in M, and the indices
%   of the maximum values in vector I.
%
% Example: a=estruct(peaks); round(max(a))==8
% Version: $Date$ $Version$ $Author$
% See also estruct, max, estruct/min

if nargin == 1
  b = [];
end
if nargin <= 2
  dim = [];
end
id=[];

% handle input estruct arrays
if numel(a) > 1 & isa(a,'estruct')
  m = zeros(size(a)); id= m;
  for index=1:numel(a)
    [m(index), id(index)] = max(a(index), b, dim);
  end
  return
end

if ~isa(a, 'estruct')
  [m,id] = max(b, a, dim);
  return
end

% return a scalar for min(a)
if isempty(b)
  m = getaxis(a, 'Signal');
  if isempty(dim), dim=1; end
  [m,id] = max(m(:),[],dim);
  return
end

% find intersection between estruct objects
cmd=a.Command;
if isa(b, 'estruct')
  [a,b] = intersect(a,b);
  m = copyobj(a);
  set(m, 'Signal', max(get(a,'Signal'), get(b,'Signal')));
  return
else
% handle estruct and scalar/vector/matrix min/max
  m = copyobj(a);
  if isempty(dim) || ~isempty(b)
    y = max(getaxis(a,'Signal'), b);
    id=[];
    set(m, 'Signal', y);
  else
    rmaxis(m); % delete all axes
    % copy all axes except the one on which operation runs
    ax_index=1;
    for index=1:ndims(a)
      if index ~= dim
        setaxis(m, ax_index, getaxis(a, num2str(index)));
        ax_index = ax_index+1;
      end
    end
    [y,id] = max(getaxis(a,'Signal'), [], dim);
    set(m, 'Signal', y, [mfilename ' of ' label(a) ]);     % Store Signal
  end
end
m.Command=cmd;
m = history(m, mfilename, a, b);


