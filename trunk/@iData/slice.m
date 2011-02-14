function slice(a)
% h = slice(s) : Plot a 3D object with slice rendering
%
%   @iData/slice function to view slices in 3D object
%     The plot is obtained with Matlab Central sliceomatic
%
% input:  s: object or array (iData)
% ex:     slice(iData(flow));
%
% Version: $Revision: 1.4 $
% See also iData, iData/plot, sliceomatic

if ndims(a) < 3 || isvector(a)
  iData_private_error(mfilename, [ 'Slice-o-matic is only available for 3D objects, but ' a.Tag ' is ' num2str(ndims(a)) '-th dimensions.' ]);
end

if prod(size(a)) > 1e6
  iData_private_warning(mfilename, [ 'Object ' a.Tag ' is too large (numel=' num2str(prod(size(a))) ...
    '.\n\tYou should rebin with e.g. a=a(1:2:end, 1:2:end, ...).' ]);
end

if exist('sliceomatic')
  x=unique(getaxis(a,2));
  y=unique(getaxis(a,1));
  z=unique(getaxis(a,3));
  c=getaxis(a,0);
  sliceomatic(c, x,y,z);
end


