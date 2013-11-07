function s = trapz(a,dim, varargin)
% s = trapz(a,dim) : computes the integral of iData objects elements along given dimension
%
%   @iData/trapz function to compute the integral of the data set along a given dimension
%     trapz(a,dim) integrates along axis of rank dim. The axis is then removed.
%       default is to use dim=1. If dim=0, integration is done on all axes and 
%       the total is returned as a scalar value. If dim='radial', the integration 
%       is done radially (R=sqrt(sum(axes^2)).
%       trapz is complementary to sum and camproj, but takes into account axes.
%     trapz(a,'radial', center) specifies the 'center' of the integration 
%       (vector of coordinates) or a single value used as center on all axes 
%       (for instance 0). All axes are assumed to be distances.
%
% input:  a: object or array (iData/array of)
%         dim: dimension to integrate (int/array of, or 'radial')
% output: s: integral of elements (iData/scalar)
% ex:     c=trapz(a);
%
% Version: $Revision$
% See also iData, iData/cumsum, iData/camproj, iData/sum

if ~isa(a, 'iData')
  iData_private_error(mfilename,[ 'syntax is ' mfilename '(iData, dim)' ]);
end

if nargin < 2, dim=1; end

if strcmp(dim, 'radial')
  s = cart2sph(a, varargin{:});
else
  s = iData_private_sumtrapzproj(a,dim, 'trapz');
end

