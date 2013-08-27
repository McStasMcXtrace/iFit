function s = kmeans(a, k)
% b = kmeans(X, k) : k-means clustering (segmentation into classes) of iData object
%
%   @iData/kmeans function to partition the object X into k classes.
%
%   b = kmeans(a,k) partitions the points in the iData object X into k clusters.
%   The resulting object Signal contains numbers from 1 to 'k' which are indices
%   of the initial data segments.
%
%   For 1D and 2D objects, the Otsu method is used, for further dimensions, the
%   legacy k-means is used (slower).
%
%   Otsu N., A threshold selection method from gray-level histogram, 
%     IEEE Trans. Syst. Man Cybern. 9:62-66;1979 
%
% input:  X: object or array (iData)
%         k: number of partitions wanted (integer, default is 2)
% output: b: object or array (iData)
% ex:     b=kmeans(a);
%
% Version: $Revision: 1035 $
% See also iData, iData/uminus, iData/abs, iData/real, iData/imag, iData/uplus

if nargin < 2
  k = 2;
end

% handle input iData arrays
if numel(a) > 1
  s = zeros(iData, numel(a),1);
  parfor index=1:numel(a)
    s(index) = feval(mfilename, a(index), k);
  end
  s = reshape(s, size(a));
  return
end

% now call Otsu of fastCMeans depending on dimensionality
X = subsref(a,struct('type','.','subs','Signal'));
if ndims(a) <= 2
  X = otsu(X, k);
else
  X = uint8(X/max(X(:))*2^8); % this is faster and requires much less memory that uint16
  X = FastCMeans(X, k);
end

% create the final object
s=copyobj(a); 
s=iData_private_history(s, mfilename, a, k);

s=set(s, 'Signal', X, 'Error', 0);
s=label(s, 0, 'Clusters/partitions');

