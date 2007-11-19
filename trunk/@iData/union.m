function [ai,bi] = union(a, b)
% [ai,bi] = union(a, b) : computes object union area and values
%
%   @iData/union function to pre-combine data sets
%   This function computes the common axes union between data sets
%   Resulting objects are returned, e.g. for performing further operations
%     ai = union(a) where 'a' is an object array computes union of all elements
%
% input:  a: object or array (iData)
%         b: object (iData)
% output: b: object or array (iData)
% ex:     b=union(a, a);
%
% See also iData, iData/setaxis, iData/getaxis, iData/intersect

if nargin == 2
  bi = interp([a b]);
  ai = bi(1);
  bi = bi(2);
  return
end

if length(a) == 1, ai=a; bi=a; return; end

% initiate new axes
for index=1:ndims(a(1))
  c_step{index} =  Inf;
  c_min{index}  =  Inf;
  c_max{index}  = -Inf;
  c_len{index}  =  0;
end

% loop on all iData to find intersection area
for index=1:length(a)
  if ndims(a(index)) ~= ndims(a(1))
    iData_private_error(mfilename, [ 'Object intersection requires same dimensionality.\n\tobject ' inputname(1) ' ' a(1).Tag ' is ' num2str(ndims(a(1))) ' but object ' a(index).Tag ' is ' num2str(ndims(a(index))) ]);
  end
  for j_ax = 1:ndims(a(index))  % for each dimension
    x = getaxis(a(index), j_ax); x=unique(x(:));    % extract axis, and remove duplicates. diff > 0
    c_step{j_ax} = min(min(diff(x)), c_step{j_ax}); % smallest step
    c_min{j_ax}  = min(min(x), c_min{j_ax});        % lowest min
    c_max{j_ax}  = max(max(x), c_max{j_ax});        % highest max
    c_len{j_ax}  = c_len{j_ax} + length(x);         % cumulated axes length
  end
end

% build new axes
for j_ax = 1:ndims(a(index))  % for each dimension
  c_len{j_ax} = c_len{j_ax}/length(a);                  % mean axis length from original data
  len         = (c_max{j_ax}-c_min{j_ax})/c_step{j_ax}; % theoretical axis length
  c_len{j_ax} = min(len, 10*c_len{j_ax});               % can not extend axes more than 10 times
  c_axis{j_ax} = linspace(c_min{j_ax}, c_max{j_ax}, c_len{j_ax});
end

% loop on all iData to interpolate
ai = a; bi=[];
for index=1:length(a)
  ai(index) = interp(a(index), c_axis);
end

