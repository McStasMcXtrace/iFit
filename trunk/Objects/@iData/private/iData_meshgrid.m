function [out, changed, flag_linspace] = iData_meshgrid(in, Signal, method)
% iData_meshgrid: checks/determine an axis system so that 
% * it is regular
% * it matches the Signal dimensions 
%
% arguments
%   in:     is a cell array of axes
%   Signal: array dimensions (as obtained from size), or array which dimensions 
%             are used for checking axes.
%   out: a cell array of axes which are regular
%   changed: true when the Signal needs to be interpolated on new axes
%   flag_linspace: array of axes that were set new vectors

if nargin < 2, Signal=[]; end
if nargin < 3, method=''; end

if isa(in, 'iData')
  if isempty(Signal), Signal = size(in); end
  out = cell(1,ndims(in));
  for index=1:length(out)
    out{index} = getaxis(in, index);  % loads object axes, or 1:end if not defined 
  end
  in = out;
end

changed = 0;

% if no Signal is defined, we will determine its size from the axes dimensions
if isempty(Signal)
  for index=1:length(in)
    if isvector(in{index});
      Signal(index) == length(in{index});
    else
      % we assume initial axes are already grid-like, and match the Signal
      Signal = size(in{index}); 
    end
  end
else
  % is Signal already a size ?
  if ~isvector(Signal) || length(Signal) ~= length(in)
    Signal = size(Signal);
  end
end

% check if axes are monotonic (that is their 'unique' values corresponds to the 
% Signal dimensions).
if length(find(Signal>1)) == 1 % vector or event type data set
  flag_linspace = ones(1, length(in));
else
  flag_linspace = zeros(1, length(in));
  for index=1:length(in)
    x = in{index}; sx = size(x); j = find(sx > 1);

    if any(sx(j) ~= Signal(j))
      % does not match Signal size (a vector or mesh with wrong size) 
      % in principle the iData internal checks prevent that but we never know...
      flag_linspace(index)=1;
    elseif ~isempty(strfind(method, 'vector')) && ~isvector(x)
      % vector requested but axis is not a vector
      flag_linspace(index)=1; 
    elseif ndims(x) == 2 && (~issorted(x, 'rows') || ~issorted(x', 'rows'))
      % axis is not sorted
      flag_linspace(index)=1; % may have good length, sorting is wrong
    end
  end
end

% determine regular vector-type axes for those flagged
out = in;
for index=1:length(in)
  x = in{index}; 
  if flag_linspace(index) == 1
    x=x(:); ux = unique(x);
    if length(ux) == Signal(index) 
      % we get a nice vector from the initial axis values
      out{index} = ux; changed=1;
    else
      % we use a new regular vector
      out{index} = linspace(min(x), max(x), Signal(index));
      changed    = 1; % new axis requires interpolation
    end
  else
    % use fake index axes just for ndgrid to work
    if numel(x) == Signal(index)
      out{index} = x;
    else
      out{index} = 1:Signal(index);
    end
  end
end

% make sure we have grid style axes
if isempty(strfind(method, 'vector'))
  [out{:}] = ndgrid(out{:});
end

