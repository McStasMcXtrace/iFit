function ax = iFunc_feval_guess_axes(model, p, varargin)

ax = varargin;
model.Dimension = abs(model.Dimension);
% complement axes if too few are given
if length(ax) < model.Dimension
  % not enough axes, but some may be given: we set them to 'empty' so that default axes are used further
  for index=(length(ax)+1):model.Dimension
    ax{index} = [];
  end
end

% default return value...
parameter_names = lower(model.Parameters);

% check axes and define missing ones
for index=1:model.Dimension
  % check for default axes to represent the model when parameters are given
  % test parameter names
  
  width    = NaN;
  position = NaN;
  for index_p=1:length(parameter_names)
    if  ~isempty(strfind(parameter_names{index_p}, 'width')) ...
      | ~isempty(strfind(parameter_names{index_p}, 'tau')) ...
      | ~isempty(strfind(parameter_names{index_p}, 'damping')) ...
      | ~isempty(strfind(parameter_names{index_p}, 'gamma'))
      if isnan(width)
        width = abs(p(index_p)); 
        % this parameter name is removed from the search for the further axes
        parameter_names{index_p} = ''; 
      end
    elseif ~isempty(strfind(parameter_names{index_p}, 'centre')) ...
      |    ~isempty(strfind(parameter_names{index_p}, 'center')) ...
      |    ~isempty(strfind(parameter_names{index_p}, 'position'))
      if isnan(position), 
        position = p(index_p);
        % this parameter name is removed from the search for the further axes
        parameter_names{index_p} = '';
      end
    end
    if ~isnan(width) && ~isnan(position)
      if isempty(ax{index})
		    % axis is not set: use default axis from parameter names and values given
		    if model.Dimension > 2, sz_max = 20; else sz_max = 50; end
		    ax{index} = linspace(position-3*width,position+3*width, sz_max+index);
		    % orient the axis along the right dimension to indicate this is not an event type
        d = ones(1,max(2,model.Dimension)); d(index) = numel(ax{index});
        ax{index} = reshape(ax{index}, d);
        width = NaN; position = NaN;
        break  % go to next axis (exit index_p loop)
      end
    end
  end % for index in parameter names
  if isempty(ax{index}) % when can not get axis, we use -5:5
    ax{index} = linspace(-5,5,50+index);
    % orient the axis along the right dimension to indicate this is not an event type
    d = ones(1,max(2,model.Dimension)); d(index) = numel(ax{index});
    ax{index} = reshape(ax{index}, d);
  end
  
end % for index in model dim

% convert axes to nD arrays for operations to take place
% check the axes and possibly use ndgrid to allow nD operations in the
% Expression/Constraint. Only for non event style axes.

% check for vectors: 0 or 1 for each axis
check_vector = cellfun(@(c)(~isempty(c) && length(c) == numel(c)), ax(1:model.Dimension));
% check number of elements: would be equal for 3D grids (xyz) or 4D (xyzt)
check_numel  = cellfun(@numel, ax(1:model.Dimension));
% check orientation (axis rank) for vectors. 0 for non vectors
check_orient = zeros(1,model.Dimension);
index        = find(check_vector);
check_orient(index) = cellfun(@(c)min(find(size(c)==numel(c))), ax(index));

if model.Dimension > 1 && all(check_vector) ...
  && all(check_numel  == check_numel(1) | check_numel == 1) ...
  && all(check_orient == check_orient(1))
  % event nD: all axes as vectors, same length, same orientation
  is_event = true;
else
  % map nD:   any other choice
  is_event = false;
end

if model.Dimension > 1
  if ~is_event && all(check_vector) && any(check_orient ~= check_orient(1))
    % when axes are all vectors, but not same orientation, we create a grid
    [ax{1:model.Dimension}] = ndgrid(ax{1:model.Dimension});
  elseif is_event
    % make sure all axes will be 'event' style, ie vectors same orientation
    sz = [];
    % first get the size of the event/cloud (first non scalar axis)
    for index=1:model.Dimension
      if ~isscalar(ax{index}), sz = size(ax{index}); break; end
    end
    if ~isempty(sz)
      % then convert all scalar stuff into same length vectors (constant)
      for index=1:model.Dimension
        if isscalar(ax{index})
          ax{index} = ax{index}*ones(sz);
        else
          ax{index} = reshape(ax{index}, sz);
        end
      end
    end
  end
end

