function slice(a)
% slice(s) : Plot a 3D object with slice rendering
%
%   @iData/slice function to view slices in 3D object
%     The plot is obtained with Matlab Central sliceomatic
%
% input:  s: object or array (iData)
% ex:     slice(iData(flow));
%
% Version: $Revision: 1.9 $
% See also iData, iData/plot, sliceomatic

if ndims(a) < 3
  iData_private_error(mfilename, [ 'Slice-o-matic is only available for 3D objects, but ' a.Tag ' is ' num2str(ndims(a)) '-th dimensions. Use plot instead.' ]);
elseif  isvector(a)
  iData_private_error(mfilename, [ 'Use hist to create a volume to display with Slice-o-matic or use plot instead.' ]);
end

if prod(size(a)) > 1e6
  iData_private_warning(mfilename, [ 'Object ' a.Tag ' is large (numel=' num2str(prod(size(a))) ...
    '.\n\tNow rebinning for display purposes with e.g. a=reducevolume(a);' ]);
  a=reducevolume(a);
end

if exist('sliceomatic')
  x=unique(getaxis(a,2));
  y=unique(getaxis(a,1));
  z=unique(getaxis(a,3));
  c=getaxis(a,0);
  sliceomatic(c, x,y,z);
  title(a.Title,'interpreter','none');
  xlabel(xlabel(a)); ylabel(ylabel(a)); zlabel(zlabel(a)); 
end


