function s = camproj(a,dim)
% s = camproj(a,dim) : computes the projection of iData objects elements
%
%   @iData/camproj function to compute the projection/sum of the elements of the data set
%     camproj(a,dim) projects along axis of rank dim. All other axes are removed.
%       If dim=0, projection is done on all axes and the total is returned as a scalar value. 
%       camproj(a,1) projects on first dimension (rows).
%       camproj is the complementary to sum.
%
% input:  a: object or array (iData/array of)
%         dim: dimension to project (int)
% output: s: projection of elements (iData/scalar)
% ex:     c=camproj(a);
%
% Version: $Revision: 1.12 $
% See also iData, iData/plus, iData/prod, iData/cumsum, iData/mean, iData/sum, iData/trapz

if ~isa(a, 'iData')
  iData_private_error(mfilename,[ 'syntax is ' mfilename '(iData, dim)' ]);
end

if nargin < 2, dim=1; end

s = iData_private_sumtrapzproj(a,dim, 'camproj');
