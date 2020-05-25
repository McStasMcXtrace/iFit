function s = cart2sph(a, center)
% CART2SPH Transform object from Cartesian to spherical coordinates.
%   CART2SPH(A) Transforms the object axes into spherical coordinates, which
%   allow radial integration aferwards.
%   The 'center' of the distribution (1st moment) is used as symmetry point 
%   for the computation of the radius. All axes are considered to be
%   distances.
%
%   CART2SPH(a,center) specifies the 'center' of the integration 
%   (vector of coordinates) or a single value used as center on all
%   axes (for instance 0).
%
% Example: a=iData(peaks); b=cart2sph(a); ndims(a) == ndims(b)
% Version: $Date$ $Version$ $Author$
% See also iData, iData/rotate, iData/sum, iData/trapz, iData/camproj

if nargin < 2, center=[]; end

% handle input iData arrays
if numel(a) > 1
  s = zeros(iData, numel(a),1);
  for index=1:numel(a)
    s(index) = feval(mfilename, a(index), center);
  end
  s = reshape(s, size(a));
  return
end

if ~isvector(a)
  a = meshgrid(a); % this way we have grid axes, vector/event objects are already OK
end

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
  warning(...
    [ mfilename ': The centroid vector is of length ' num2str(length(center)) ' but the object requires ' num2str(ndims(a)) ' values (dimensionality).' ]);
  return;
end

s = copyobj(a);
cmd= a.Command;

% we interpolate the object on a grid so that all axes have the same size
if ndims(a) == 2
  y = private_cleannaninf(getaxis(a, 1)) - center(1);
  x = private_cleannaninf(getaxis(a, 2)) - center(2);
  [theta, rho] = cart2pol(x,y);
  phi = [];
  clear y
  x = getaxis(a, 0);
elseif ndims(a) == 3
  y = private_cleannaninf(getaxis(a, 1)) - center(1);
  x = private_cleannaninf(getaxis(a, 2)) - center(2);
  z = private_cleannaninf(getaxis(a, 3)) - center(3);
  [theta, phi, rho] = cart2sph(x,y,z);
  clear y z
  x = getaxis(a, 0);
else
  rho = zeros(size(s));
  % we extract Signal and all axes, except 'dim'
  for index=1:ndims(a)
    % then compute the sqrt(sum(axes.*axes))
    x            = private_cleannaninf(getaxis(a, index)) - center(index);
    rho          = rho + x.*x;
    clear x
  end
  theta=[]; phi=[];
  rho = sqrt(rho);
  x = getaxis(a, 0);
end

% create the output object

% Store Signal and Monitor
s = setalias(s, 'Signal',  x);

s = rmaxis(s);  % remove all axes, will be rebuilt after operation

s = setaxis( s, 1, rho);
s = label(   s, 1, 'Radius');

if ~isempty(theta)
  s = setaxis( s, 2, theta*180/pi);
  s = label(   s, 2, 'Theta Azimutal [deg]');
end
if ~isempty(theta)
  s = setaxis( s, 3, phi*180/pi);
  s =  label(  s, 3, 'Phi Elevation [deg]');
end
s = setalias(s, 'Center', center);
s = history(s, mfilename, a, center);

