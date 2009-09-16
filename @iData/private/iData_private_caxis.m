function c_axis=iData_private_caxis(a)
% compute common axis for union

c_axis=[];
if isempty(a), return; end

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
    y = min(min(diff(x)), c_step{j_ax}); % smallest step
    if ~isempty(y), c_step{j_ax}=y; end    
    c_min{j_ax}  = min(min(x), c_min{j_ax});        % lowest min
    c_max{j_ax}  = max(max(x), c_max{j_ax});        % highest max
    c_len{j_ax}  = c_len{j_ax} + length(x);         % cumulated axes length
  end
end

% build new axes
for j_ax = 1:ndims(a(index))  % for each dimension
  c_len{j_ax} = c_len{j_ax}/length(a);                  % mean axis length from original data
  len         = (c_max{j_ax}-c_min{j_ax})/c_step{j_ax}; % theoretical axis length
  c_len{j_ax} = min(len+1, 10*c_len{j_ax});             % can not extend axes more than 10 times
  c_axis{j_ax} = linspace(c_min{j_ax}, c_max{j_ax}, c_len{j_ax});
end

