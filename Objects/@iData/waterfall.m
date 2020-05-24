function h = waterfall(a, option)
% WATERFALL Waterfall plot (2D/3D).
%   H = WATERFALL(S,OPTION) plot a 2D/3D object as waterfall (side by side 
%   lines). The plot graphic object handle H is returned.
%
% Example: h=waterfall(iData(flow)); tf=ishandle(h); delete(gcf); tf
% Version: $Date$ $Version$ $Author$
% See also iData, iData/plot

if nargin ==1
	option='';
end
h = plot(a, [ 'waterfall ' option ]);

