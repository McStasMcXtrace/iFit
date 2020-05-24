function b = interp(self, varargin)
% INTERPN N-D object interpolation (table lookup).
%    VI = INTERPN(V, Y1,Y2,Y3,...) interpolates object V onto new axes Y1,Y2,Y3, 
%    etc to find the underlying object VI values.
%    INTERPN works for all object dimensionalities. The multidimensional interpolation
%    is based on a Delauney tessellation using the Computational Geometry 
%    Algorithms Library, CGAL. Extrapolated data is set to NaN for the Signal, 
%    Error and Monitor. Also, in some cases, the triangulation interpolant creates
%    fake 'flat' area, especially when the axes area is concave. We then recommend
%    you try the HIST method.
%    For Event data sets, we recommend to use the HIST method which is much faster.
%
%    VI = INTERPN(V, {Y1,Y2,Y3,...}) is similar to the previous syntax.
%
%    VI = INTERPN(V, A) where A is an iData object interpolates V onto A axes,
%    i.e. Y1=A{1}, Y2=A{2}, ... with the usual iData axis notation.
%
%    VI = INTERPN(...,METHOD) specifies alternate methods.  The default
%    is linear interpolation.  Available methods are:
% 
%      'nearest' - nearest neighbor interpolation
%      'linear'  - linear interpolation (default)
%      'spline'  - spline interpolation
%      'cubic'   - cubic interpolation as long as the data is uniformly
%                  spaced, otherwise the same as 'spline'
%
%    VI = INTERPN(..., 'grid') uses meshgrid/ndgrid to determine new axes as arrays
%
% Example: a=iData(peaks); b=interp(a, 'grid'); isequal(a,b)
% Version: $Date$ $Version$ $Author$
% See also iData, interp1, interpn, ndgrid, iData/setaxis, iData/getaxis,
%          iData/hist, iData/resize, iData/reshape, iData/fill

b = interpn(self, varargin{:});
