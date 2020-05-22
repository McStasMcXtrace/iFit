function [b, j] = jacobian(a, varargin)
% JACOBIAN Compute the Jacobian.
%   [B,J]=JACOBIAN(S, D) transforms 'S' on the 'D' axes and computes the Jacobian 
%   from the object current axes to new axes. This is to be used for object 
%   coordinate/variable changes from the given object axes to new axes. The 
%   'from' and 'to' axes should have the same dimensions.
%   The transformed object is returned as B, as well as the Jacobian J which was  
%   used for the transformation.
%   The Jacobian is not a mere axis assignement, nor an interpolation.
%   The Jacobian transform retains the full integral of the initial object, that
%   is:
%      trapz(s,0) == trapz(jacobian(s, ...),0)
%
%   A direct axes assignation is perfomed by using 'setaxis' (which does not affect
%   the Signal, and retains its sum, but not its integral). 
%   A direct axes rebinning can be performed by using 'interp' (which affects the
%   Signal, its sum, but usually retains its integral).
%
%   B=JACOBIAN(S, X1,X2, ... Xn) transforms S onto 'X1...Xn' axes given as 
%   vectors or matrices as obtained from ndgrid, up tu NDIMS(S).
%
%   B=JACOBIAN(S, {X1,X2, ... Xn}) is the same as above.
%
%   B=JACOBIAN(S, ..., 'xlab','ylab',...) specifies the new axes labels.
%
%   B=JACOBIAN(S, ..., {'xlab','ylab',...}) is the same as above.
%
% Example: a=estruct(peaks); x=linspace(1,2,size(a,1)); ...
%          g=jacobian(a, x, [],'half X'); trapz(a,0)==trapz(g,0)
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/del2, diff, estruct/gradient, estruct/interp, 
%   estruct/setaxis, gradient

% TODO: this implementation computes the determinent assuming axes are quasi-perpendicular
% then det = product(trace) as done in the function gradient_axis.

if nargin <= 1,
  b = a; j=[]; return
end

% handle input estruct arrays
if numel(a) > 1
  b =zeros(estruct, size(a));
  for index=1:numel(a)
    b(index) = feval(mfilename,a(index), varargin{:});
  end
  return
end

% 'from' axes
i_axes = cell(1,ndims(a)); i_labels=i_axes;
for index=1:ndims(a)
  [i_axes{index}, i_labels{index}] = getaxis(a, index);  % loads object axes, or 1:end if not defined 
end

% 'to' axes
f_axes           = i_axes;
f_labels         = i_labels;

% parse varargin to overload defaults and set manually the axes
axis_arg_index   = 0;
label_arg_index  = 0;
for index=1:length(varargin)
  c = varargin{index};
  if isa(varargin{index}, 'estruct')  % get 'to' axis from other estruct object
    if length(c) > 1
      warning([ mfilename, ...
          ': Can not apply Jacobian onto all axes of input argument ' num2str(index) ...
           ' which is an estruct array of ' num2str(numel(c)) ...
           ' elements. Using first element only.']);
      c = c(1);
    end
    for j1 = 1:ndims(c)
      axis_arg_index = axis_arg_index+1;
      [f_axes{axis_arg_index}, lab] = getaxis(c, j1);
      if ~isempty(lab) && axis_arg_index < length(i_labels) && isempty(i_labels{axis_arg_index})
        i_labels{axis_arg_index} = lab;
      end
    end
  elseif isnumeric(c) & ~isempty(c)   % get 'to' axis from  vector/matrix
    axis_arg_index = axis_arg_index+1;
    f_axes{axis_arg_index} = c;
  elseif ischar(c)                    % get final label
    label_arg_index = label_arg_index+1;
    if ~isempty(c) f_labels{label_arg_index} = c; end
  elseif iscell(c) && ischar(c{1})    % get final label (as cellstr)
    for j1 = 1:length(c(:))
      label_arg_index = label_arg_index+1;
      if ~isempty(c{j1}) f_labels{label_arg_index} = c{j1}; end
    end
  elseif iscell(c)                    % get 'to' axis from  cell(vector/matrix)
    for j1 = 1:length(c(:))
      axis_arg_index = axis_arg_index+1;
      if ~isempty(c{j1}) f_axes{axis_arg_index} = c{j1}; end
    end
  elseif ~isempty(c)
    warning([ mfilename, ...
        ': Input argument ' num2str(index) ' of class ' class(c) ' size [' ...
         num2str(size(c)) '] is not supported. Ignoring.']);
  end
end
clear varargin

% tests axes
if ndims(a) > 1
  for index=1:ndims(a)
    % test for the target axes in case they are given as scalars (axes spacing)
    if isscalar(f_axes{index}) && ~isscalar(i_axes{index})
      x = i_axes{index}; x=unique(x); % also makes it a vector
      f_axes{index} = min(x):f_axes{index}:max(x);
    end
    myisvector = @(c)max(size(c)) == numel(c);
    if myisvector(f_axes{index}) % vectors should be oriented the right way
      d=ones(1, ndims(a));
      d(index) = length(f_axes{index});
      f_axes{index} = reshape(f_axes{index}, d);
    end
    try
      % axes must have same size as the initial ones
      if ~all(size(i_axes{index}) == size(f_axes{index}))
        warning([ mfilename,...
            ': Jacobian can only be computed if axes have same size. Object axis{%i}=[' ...
             num2str(size(i_axes{index})) '], new axis=[' size(f_axes{index}) '].' ]);
        
      end 
    end
  end
end

% compute the Jacobian product for each initial axis
gi = gradient_axis(i_axes);

% compute the Jacobian product for each final axis
gf = gradient_axis(f_axes);

% now divide the gradients to form the Jacobian
cmd=a.Command;
b = copyobj(a);
j = genop(@rdivide, gf, gi);
j = abs(j);
clear gi gf
s = get(b,'Signal')./j;
e = get(b,'Error') ./j;
if nargout < 2
  clear j
elseif all(j(:) == j(1))
  j = j(1);
end

% assemble the resulting object
[dummy, sl] = getaxis(b, '0');
b = set(b, 'Signal', s, 'Error', abs(e));
label(b, 0, [  mfilename '(' sl ')' ]);

% redefine axes
for index=1:length(f_axes)
  if isequal(i_axes{index},f_axes{index})
    continue; %axes are the same...
  else
    ax=[ 'Axis_' num2str(index) ];
    setalias(b, ax, f_axes{index});
    setaxis(b, index, ax);
  end
  label(b, ax, f_labels{index});
end

b.Command=cmd;
b = history(b, mfilename, a);  

end
% ------------------------------------------------------------------------------
function gp=gradient_axis(ax)
  gp = [];
  myisvector = @(c)length(c) == numel(c);
  for index=1:length(ax)
    s = ax{index};
    grad_length = ndims(s);
    if myisvector(s), grad_length=1; end
    g = cell(1,grad_length);
    [g{:}]=gradient(s); % compute the gradient there
    gs = [];
    % now computes the product(g) for a given gradient(axis)
    for dim=1:grad_length
      if ~isempty(g{dim}) && any(g{dim})
        if isempty(gs), gs = g{dim};
        else            gs = genop(@times, gs, g{dim});
        end
      end
    end % for dim
    % orient vector axes along the right dimension
    if myisvector(gs)        
      v        = ones(1,length(size(s)));
      v(index) = length(gs);
      gs       = reshape(gs, v);
    end
    % gradient product for all dimensions
    if isempty(gp), gp = gs;
    else            gp = genop(@times, gp, gs);
    end
  end % for index
end


