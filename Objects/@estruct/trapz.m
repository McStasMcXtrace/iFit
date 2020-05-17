function [s,sigma] = trapz(a,dim, varargin)
% TRAPZ  Trapezoidal numerical integration.
%   S = TRAPZ(A,DIM) integrates along axis of rank DIM. The axis is then removed.
%   Default is to use DIM=1. 
%   TRAPZ is complementary to SUM and CAMPROJ, but takes into account axes values.
%
%   [S,err] = TRAPZ(A,0) integrates along all axes and the total is returned as a 
%   scalar value. In this case, a second output argument holds the error bar.
%
%   S = TRAPZ(A,'radial') integrates radially (R=sqrt(sum(axes^2)).
%
%   TRAPZ(A,'radial', center) specifies the 'center' of the integration 
%   (vector of coordinates) or a single value used as center on all axes 
%   (for instance 0). All axes are assumed to be distances.
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/cumsum, estruct/camproj, estruct/sum, estruct/cart2sph

if nargin < 2, dim=1; end

if strcmp(dim, 'radial')
  sigma = [];
  s = cart2sph(a, varargin{:});
else
  [s,sigma] = private_sumtrapzproj(a,dim, 'trapz');
end

