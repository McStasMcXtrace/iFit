function h=plot(a, varargin)
% PLOT Display a plot of an object.
%   H = PLOT(S) plot N-dimensional object, i.e. the signal of the object as a
%   function of the defined axes. The plot is an errorbar plot for 1D data y=f(x),
%   a surface mesh for 2D z=f(x,y), and an isosurface for 3D c=f(x,y,z) data.
%   Data specified as series of coordinate points are plotted using a plot3-type
%   rendering. Object with dimensions above 3 are projected into a 3D space.
%   The resulting graphic handle is returned into H.
%
%   H = PLOT(S1, S2, ...) plot multiple objects onto the same coordinate frame.
%   This corresponds to an overlay-type plot. This syntax is equivalent to
%   H = PLOT([ S1 S2 S3 ...]).
%
%   H = PLOT(..., METHOD) specifies a rendering option.
%   For 1D plots y=f(x), METHOD is a string to specify color/symbol, such
%     as in the usual PLOT function, e.g. METHOD='r-'. The 'hide_errorbars' method
%     can be added to the line/color specification in order not to plot error bars.
%     To plot a set of 1D objects side by side, specify a 2D plot option such as
%     'surf' or 'plot3'. Alternatively you may use CAT to concatenate 1D objects
%     into a 2D one.
%   For 2D plots z=f(x,y), METHOD is a string which may contain:
%     surf, mesh, contour, contour3, surfc, surfl, contourf
%     plot3, scatter3 (colored points), stem3, pcolor, waterfall
%   For 3D plots c=f(x,y,z), METHOD is a string which may contain:
%     plot3 (semi transparent volume),
%     scatter3 (colored points, supports 'scatter3 filled' and 'scatter3 bubble')
%     waterfall (supports 'waterfall x' 'y' and 'z'), contour (set of contour plots)
%     surf, surf median, surf mean, surf half (isosurface)
%
%   METHOD may also contain global options for 2D and 3D plots: 
%     flat, interp, faceted (for shading), view2, view3
%     transparent, light, clabel, colorbar, shifted (overlayed 2D)
%
%   METHOD may also contain global options for all plots: 
%     axis tight, axis auto, hide_axes (compact layout)
%     painters (bitmap drawing), zbuffer (vectorial drawing)
%     opengl (faster for large data sets)
%     whole or full (do not reduce large object size for plotting)
%     figure (open a new figure window)
%     replace (replace existing plots for same objects)
%     legend (plot legend for objects)
%
%   H = PLOT(..., METHOD, ...) additional arguments are sent to the plot command.
%
%   When a XtickLabel alias/property exists in the object, the Xtick labels are set.
%   The SLICE method opens the interactive sliceomatic 3D viewer.
%
%   As mentioned in the estruct axis definition (see estruct/setaxis), the 'X' axis
%     on plots refers to the 2nd dimension (along columns), whereas the 'Y' axis 
%     refers to the first dimension (along rows).
%
%  Type <a href="matlab:doc(estruct,'Plot')">doc(estruct,'Plot')</a> to access the iFit/Plot Documentation.
%
% Contributed code (Matlab Central): 
%   fscatter3:   Felix Morsdorf, Jan 2003
%   vol3d:       Joe Conti, 2004
%   sliceomatic: Eric Ludlam 2001-2008
%
% Version: $Date$ $Version$ $Author$
% See also estruct, interp1, interpn, ndgrid, plot, estruct/setaxis, estruct/getaxis
%          estruct/xlabel, estruct/ylabel, estruct/zlabel, estruct/clabel, estruct/title
%          shading, lighting, surf, estruct/slice, estruct/contour3, estruct/contourf
%          estruct/contour, estruct/surf, estruct/slice, estruct/plot3, estruct/surfl, 
%          estruct/surfc, estruct/mesh

h      = [];
funcs  = []; % additional iFunc objects to plot afterwards...
method = '';
args   = {};
vargs  = {};

if ~isempty(get(0,'CurrentFigure'))
  fig = get(0,'CurrentFigure');
  if ~strcmp(get(fig, 'HandleVisibility'),'on') % make sure we do not overwrite in a protected figure
    method = 'figure ';
  end
  ih     = ishold;
else ih=0; end

% analyze input arguments
if nargin == 1, 
  h = private_plot(a, method);
elseif length(varargin) == 1
  if ischar(varargin{1})
    method=[ method varargin{1} ];
    args = method;
  elseif isa(varargin{1},class(a))
    b = varargin{1};
    a = [ a(:) ; b(:) ];
  elseif isa(varargin{1},'iFunc')
    funcs = varargin{1};
  end
  h = private_plot(a, method);
else % multiple plot/methods to render
  % first extract non char/object/iFunc
  for index=1:length(varargin)
    if ~ischar(varargin{index}) && ~isa(varargin{index},class(a)) && ~isa(varargin{index},'iFunc')
      args{end+1}    = varargin{index};
      varargin{index}='';
    else
      vargs{end+1} = varargin{index};
    end
  end
  varargin =vargs; % now only contains char/object/iFunc. 'args' contains the rest
  % split varargin looking for chars
  index=1;  
  while index <= length(varargin)  % parse input arguments and split with char/methods calls
    if ~isempty(varargin{index})
        if ischar(varargin{index})
          method1 = [ method varargin{index} ]; % plot stored objects with current method
          hh = private_plot(a, method1, args{:});
          h =[ h(:) ;  hh(:) ];
          a = [];
          hold on
        elseif isa(varargin{index},class(a)) % store some objects until we plot them
          b = varargin{index};
          a = [ a(:) ; b(:) ];
        elseif isa(varargin{index},'iFunc') % store some objects until we plot them (at the end)
          funcs = [ funcs ; varargin{index} ];
        end
    end
    index=index+1;
  end
  if numel(a) > 0 % if we have some objects left, we plot them with default method
    hh = private_plot(a, '', args{:})
    h =[ h(:) ; hh(:)  ];
  end
  return
end
clear varargin

% handle legend stuff
if iscell(h)
  try
    index = strfind(args, 'legend');
    if iscell(index), index = find(~cellfun(@isempty, index)); end
  catch
    index = [];
  end
  if ~isempty(index)
    legend([ h{:} ]);
  end
end

% handle iFunc objects
if ~isempty(funcs)
  hold on
  axis(axis); % fix plot limits
  hline = plot(funcs);
  set(findobj(hline,'Type','Line'),'LineStyle','--');
  h = [ h(:) ;  hline(:) ];
end

if ~ih, hold off; end
