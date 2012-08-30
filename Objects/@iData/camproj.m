function s = camproj(a,dim, radial)
% s = camproj(a,dim) : projection/radial integration of iData objects elements
%
%   @iData/camproj function to compute the projection/sum of the elements of the data set
%     camproj(a,dim) projects along axis of rank dim. All other axes are removed.
%       If dim=0, projection is done on all axes and the total is returned as a scalar value. 
%       camproj(a,1) projects on first dimension (rows).
%       camproj is the complementary to sum.
%
%     camproj(a,dim, center) : computes the radial integration of axis 'dim' along 
%       other axes the 'center' argument can be a vector specifying center of the 
%       integration or entered as the string 'auto' to determine automatically the 
%       center from the distributions 1st moment (std, mean), or a single value used 
%       as center on all axes.
%       The 'dim' input argument is usually 0, in which case 'radial' is [y0, x0, z0, ...]
%       When 'dim' is not 0, radial is [ Signal center, axis centers except the 'dim one, ... ]
%
% input:  a:     object or array (iData/array of)
%         dim:   dimension rank to project to (int)
%         center:'auto' or a vector which length is the object dimensionality 
% output: s: projection of elements (iData 1D/scalar)
% ex:     c=camproj(a);
%
% Version: $Revision: 1.12 $
% See also iData, iData/rotate, iData/sum, iData/trapz

if ~isa(a, 'iData')
  iData_private_error(mfilename,[ 'syntax is ' mfilename '(iData, dim)' ]);
end

if nargin < 2, dim=1; end

if nargin <= 2
  s = iData_private_sumtrapzproj(a,dim, 'camproj');
else
  % radial integration (works for surfaces, spheres, ...
  if ndims(a) < 2, s=a; return; end
  
  % handle input iData arrays
  if numel(a) > 1
    s = [];
    for index=1:numel(a)
      s = [ s feval(op, a(index), dim, 'radial') ];
    end
    s = reshape(s, size(a));
    return
  end
  
  % handle center of the integration
  if ischar(radial) || isempty(radial)
    % use 1st moment for each integration axis (automatic)
    radial=zeros(1,ndims(a));
    radial_index = 1;
    for index=0:ndims(a)
      if index == dim, continue; end
      [dummy, radial(radial_index)] = std(a, index);
      radial_index = radial_index+1;
    end
  end
  
  if isnumeric(radial) && length(radial) == 1
    radial = radial*ones(1,ndims(a));
  end
  
  if length(radial) < ndims(a)
    iData_private_warning(mfilename, ...
      [ 'The centroid vector is of length ' num2str(length(radial)) ' but the object requires ' num2str(ndims(a)) ' values (dimensionality).' ]);
    return;
  end
  
  s = copyobj(a);
  cmd= a.Command;

  % we interpolate the object on a grid so that all axes have the same size
  a = interp(a, 'grid');
  if ndims(a) == 2
    index=0:2; index(dim+1) = [];
    x = iData_private_cleannaninf(getaxis(a, index(1))) - radial(1);
    y = iData_private_cleannaninf(getaxis(a, index(2))) - radial(2);
    rho = hypot(x,y); % faster and more accurate
    clear x y
  else
    rho = zeros(size(s));
    radial_index=1;
    % we extract Signal and all axes, except 'dim'
    for index=0:ndims(a)
      % then compute the sqrt(sum(axes.*axes))
      if index == dim, continue; end
      x            = iData_private_cleannaninf(getaxis(a, index)) - radial(radial_index);
      radial_index = radial_index+1;
      rho          = rho + x.*x;
      clear x
    end
    rho = sqrt(rho);
  end
  
  % create the output object
  x = getaxis(a, dim);
  [rho, index] = sort(rho(:));
  x=x(:); x = x(index);
  % Store Signal and Monitor
  s = setalias(s, ...
    'Signal',  x(:), [ 'radial integration of ' label(a, dim) ' axis ' num2str(dim) ]);
  s = setalias(s, ...
    'Monitor', get(a,'Monitor'));
  s = rmaxis(s);  % remove all axes, will be rebuilt after operation
  if dim ~= 0
    s = set(s, 'Error', 0);
  end
  s = setaxis(s, 1, rho);
  s = xlabel(s, 'Radius');
  s = setalias(s, 'Center', radial);
end
