function h =  pcolor(a, option)
% PCOLOR Pseudocolor (checkerboard) plot.
%   PCOLOR(S) is a pseudocolor or "checkerboard" plot of 2D object S.
%   The values of the elements of C specify the color in each
%   cell of the plot. S is displayed as a flat matrix image.
%
%   H =  PCOLOR(S,OPTION) specifies plot options for 2D and 3D plots: 
%                 flat, interp, faceted (for shading)
%                 transparent, light, clabel
%                 axis tight, axis auto, view2, view3
%                 painters (bitmap drawing), zbuffer (vectorial drawing)
%
% Example: h=pcolor(estruct(peaks)); tf=ishandle(h); delete(gcf); tf
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plot

if nargin ==1
	option='';
end
h = plot(a, [ ' pcolor ' option ]);



