function s = camproj(a,dim, center)
%  CAMPROJ Projection/radial integration of objects elements.
%    CAMPROJ(a,dim) projects along axis of rank dim. All other axes are removed.
%    CAMPROJ is the complementary to SUM.
%
%    CAMPROJ(a,0) projection is done on all axes and the total is returned as a 
%    scalar value. CAMPROJ(a,1) projects on first dimension (rows).
%
%    CAMPROJ(a,'radial') computes the radial integration (R=sqrt(sum(axes^2)).
%    the 'center' of the distribution (1st moment) is used as symmetry point 
%    for the computation of the radius.
%
%    CAMPROJ(a,'radial', center) specifies the 'center' of the integration 
%    (vector of coordinates) or a single value used as center on all axes 
%    (for instance 0). All axes are assumed to be distances.
%    The radial distribution can then be transformed into an histogram with
%    e.g. HIST(CAMPROJ(a), 100);
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/rotate, estruct/sum, estruct/trapz, estruct/cart2sph, estruct/hist

if nargin < 2, dim=1; end

if nargin <= 2 && isnumeric(dim)
  s = private_sumtrapzproj(a,dim, 'camproj');
else
  % radial integration (works for surfaces, spheres, ...)
  dim = 'radial';
  if ndims(a) < 2, s=a; return; end
  if nargin < 3, center=[]; end
  
  % handle input estruct arrays
  if numel(a) > 1
    s = zeros(estruct, numel(a),1);
    for index=1:numel(a)
      s(index) = feval(mfilename, a(index), dim, center);
    end
    s = reshape(s, size(a));
    return
  end
  
  a = meshgrid(a);
  % handle center of the integration
  if ischar(center) || isempty(center)
    % use 1st moment for each integration axis (automatic)
    center=zeros(1,ndims(a));
    for index=1:ndims(a)
      [dummy, center(index)] = std(a, index);
    end
  end
  
  if isnumeric(center) && length(center) == 1
    center = center*ones(1,ndims(a));
  end
  
  if length(center) < ndims(a)
    warning( ...
      [ mfilename ': The centroid vector is of length ' num2str(length(center)) ' but the object requires ' num2str(ndims(a)) ' values (dimensionality).' ]);
    return;
  end
  
  s = copyobj(a);
  cmd= a.Command;

  % we interpolate the object on a grid so that all axes have the same size
  if ndims(a) == 2
    x = private_cleannaninf(getaxis(a, 1)) - center(1);
    y = private_cleannaninf(getaxis(a, 2)) - center(2);
    rho = hypot(x,y); % faster and more accurate
    clear x y
  else
    rho = zeros(size(s));
    % we extract Signal and all axes, except 'dim'
    for index=1:ndims(a)
      % then compute the sqrt(sum(axes.*axes))
      x            = private_cleannaninf(getaxis(a, index)) - center(index);
      rho          = rho + x.*x;
      clear x
    end
    rho = sqrt(rho);
  end
  
  % create the output object
  x = getaxis(a, 0);
  [rho, index] = sort(rho(:));
  x=x(:); x = x(index);
  % Store Signal and Monitor
  s = setalias(s, ...
    'Signal',  x(:), [ 'radial integration of ' label(a, 0) ]);

  s = rmaxis(s);  % remove all axes, will be rebuilt after operation

  s = setaxis( s, 1, rho);
  s = xlabel(  s, 'Radius');
  s = setalias(s, 'Center', center);
  s = history(s, mfilename, a, dim, center);
end
