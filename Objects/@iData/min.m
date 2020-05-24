function [m,id] = min(a,b, dim)
% MIN    Smallest component.
%   MIN(A) returns a single value as the minimum value of the object signal.
%
%   MIN(A,B) returns an object which signal is the highest of A and B.
%
%   MIN(A,[], DIM) returns minimum value along dimension 'DIM'
%
%   [M,I] = MIN(A,...) returns the minimum value/object in M, and the indices
%   of the minimum values in vector I.
%
% Example: a=iData(peaks); round(min(a))==-7
% Version: $Date$ $Version$ $Author$
% See also iData, max, iData/min

if nargin == 1
  b = [];
end
if nargin <= 2
  dim = [];
end
id=[];

% handle input iData arrays
if numel(a) > 1 & isa(a,'iData')
  m = zeros(size(a)); id= m;
  for index=1:numel(a)
    [m(index), id(index)] = min(a(index), b, dim);
  end
  return
end

if ~isa(a, 'iData')
  [m,id] = min(b, a, dim);
  return
end

% return a scalar for min(a)
if isempty(b)
  m = getaxis(a, 'Signal');
  if isempty(dim), dim=1; end
  [m,id] = min(m(:),[],dim);
  return
end

% find intersection between iData objects
cmd=a.Command;
if isa(b, 'iData')
  [a,b] = intersect(a,b);
  m = copyobj(a);
  set(m, 'Signal', min(get(a,'Signal'), get(b,'Signal')));
  return
else
% handle iData and scalar/vector/matrix min/min
  m = copyobj(a);
  if isempty(dim) || ~isempty(b)
    y = min(getaxis(a,'Signal'), b);
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
    [y,id] = min(getaxis(a,'Signal'), [], dim);
    set(m, 'Signal', y, [mfilename ' of ' label(a) ]);     % Store Signal
  end
end
m.Command=cmd;
m = history(m, mfilename, a, b);


