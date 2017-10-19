function [b, mask, f] = imroi(a, options)
% imroi: define a region of interest on a data set
% 
% b = imroi(a)
%  This function allows to select a region-of-interest (ROI) over an existing 
%    data set, which defines an area where data points are selected. The selected 
%    data set is returned, with NaN's elsewhere.
%
% The mouse is used to define points/lines in the data set. These vertices are used
%   as a polygon shape which intersection with the data set defined the selection.
% When used with surfaces and volumes, you may orient the view (using the Rotate 
%   icon from the toolbar or Tools menu) prior the ROI selection.
%
% [b, mask] = imroi(a)
%  Same as above, but returns the mask data set, which contains 0 and 1.
%
% [b, mask, f] = imroi(a)
%  Same as above, but returns the selection figure handle.
%
% [b, mask] = imroi(a, options)
%  The options are used to customize the plot rendering, see iData/plot.
%
% Interaction:
%     mouse left-click           add a point/line
%     right-click or BACKSPACE   removes last point
%     DEL or "c"                 removes all points
%     return or "q" or middle-click: terminate input
%     ESC                        abort
%     h                          display a help dialogue
%
% Example:
%  a = iData(peaks); b=imroi(a);
%
%  input:
%     a: a data set (iData)
%
%  output:
%     b:    the data set with selected area (iData)
%     mask: the mask used for the selection, with 0 and 1 as Signal (iData)
%     f:    figure used to define the selection
%
% Version: $Date$
% See also iData, iData/plot, iData/edit
  
  if numel(a) > 1
    warning([ mfilename ': can only be used with a single data set. Using first.' ])
    a = a(1);
  end
  
  if nargin < 2, options = []; end
  
  % first we plot the object
  f = figure;
  h = plot(a, options);
  
  % we make sure the object is a meshgrid
  a = meshgrid(a);
  b = []; mask = [];
  
  % extract the axes and signal
  ax_signal = {};
  for index=[ 1:ndims(a) 0 ]
    ax_signal{end+1} = getaxis(a, index);
  end
  
  % then we request the selection (private imroi1 below)
  [select, xselect, yselect, zselect, cselect, h] = imroi1(h, ax_signal{:});
  if isempty(select), return; end
  
  % then re-build a new object with the selection
  [~, sl] = getaxis(a, '0');  % signal definition/label
  b     = copyobj(a);
  s_a   = ax_signal{end};
  clear ax_signal
  
  s_b   = s_a*nan;
  s_b(select) = s_a(select);
  b = setalias(b, 'Signal', s_b, [  'imroi(' sl ')' ]);
  b = setalias(b, 'Error',   0);
  b = setalias(b, 'Monitor', 0);
  
  % then build the mask
  if nargout > 1
    mask = copyobj(b);
    s_b = s_a*0;
    s_b(select) = 1;
    mask = setalias(mask, 'Signal', s_b, [  'imroi(' sl ') mask' ]);
  end
    
end % imroi

function [select, xselect, yselect, zselect, cselect, h] = imroi1(h, xdata, ydata, zdata, cdata)
% imroi: define a region of interest on an axis
% 
%  This function allows to select over an existing plot a set of points/lines 
%  which define an area where data points are selected.
%  You may specify which data set to use, or the active plot will be used.
%
%  input:
%   mouse left-click                add a point/line
%   right-click or DEL or BACKSPACE removes last point
%   return or 'q' or middle-click   terminate input
%   ESC                             abort
%
%  return:
%   the selected area/data set
  
  if nargin == 0, h = []; end
  if nargin < 2, xdata = []; end
  if nargin < 2, ydata = []; end
  if nargin < 4, zdata = []; end
  if nargin < 5, cdata = []; end
  xselect = []; yselect = []; zselect = []; cselect = []; select = [];
  
  if isempty(h), h = gca; end
  
  % get_children: get all child objects (in case we have an hggroup or patch)
  h = get_children(h);  % private below
  
  % get the data from the plot if not given as input or left empty
  if isempty(xdata)
    try
      xdata = get(h, 'XData');
    catch
      error([ mfilename ': the plot must have attached xdata,ydata and optionally zdata, cdata. This is not the case for active object.' ])
    end
  end
  
  if isempty(ydata)
    ydata = get(h, 'YData');
  end
  if isempty(zdata)
    try
      zdata = get(h, 'ZData');
    end
  end
  if isempty(cdata)
    try
      cdata = get(h, 'CData');
    end
  end
  if ~iscell(xdata)
    xdata = { xdata }; ydata = { ydata }; 
    zdata = { zdata }; cdata = { cdata }; 
  end
  
  % acquire points over the axes (polygon as a 'tube' seen from current orientation)
  ROIs = selectpoints; % get the selection from mouse clicks on axis
  if isempty(ROIs) || ~iscell(ROIs), return; end
  
  % loop on all available data plots with data and catenate selected points
  for index=1:numel(xdata)
  
    % get the current data set
    xd1 = xdata{index};  yd1 = ydata{index};
    zd1 = zdata{index};  cd1 = cdata{index};
    in  = [];

    for s_index = 1:numel(ROIs)  % loop on polygons
      s = ROIs{s_index};
      if ~isfield(s, 'XV')
        disp(s)
        continue
      end

      % inhull: get data points enclosed
      warning('off','inhull:degeneracy')
      % inhull below uses n x p points in p dimensions
      if ~isempty(zd1) && ~isempty(s.ZV)
        testpts = [   xd1(:)   yd1(:)   zd1(:) ];
        xyz     = [ s.XV(:)  s.YV(:)  s.ZV(:) ];
      else
        testpts = [   xd1(:)   yd1(:) ];
        xyz     = [ s.XV(:)  s.YV(:) ];
      end
      try
        in1 = inhull(testpts, xyz);
      catch
        in1 = [];
      end
      if isempty(in1), continue; end
      if isempty(in), in = in1;
      else            in = in | in1; end
    end % for index ROI
    
    % we highlight the selected region by overlaying a plot
    xsel = xd1(in); ysel = yd1(in);
    hold on
    if ~isempty(zd1) && ~isempty(s.ZV)
      zsel = zd1(in);
      plot3(xsel, ysel, zsel,'ro');
    else
      zsel = [];
      plot(xsel, ysel, 'ro');
    end
    csel = [];
    if ~isempty(cd1)
      try
        csel = cd1(in);
      end
    end
      
    % in case we have a set of objects (as a cell), we catenate all points
    if numel(xdata) > 1
      % concatenate from different objects. Can not retain the shape of the data.
      xselect = [ xselect ; xsel(:) ];
      yselect = [ yselect ; ysel(:) ];
      zselect = [ zselect ; zsel(:) ];
      select  = [ select  ; in(:)   ];
    else
      % single object. We retain the shape
      xselect = xsel; yselect = ysel; zselect = zsel; cselect = csel;
      select = in;
    end
  
  end % for index xdata
end
  
% ------------------------------------------------------------------------------  
%                              private sub-functions 
% ------------------------------------------------------------------------------

function s = selectpoints()
  % we record selections until 'end', 'abort' or 'rotate'
  % when in 'rotate' mode, we activate rotation, and then go back to selection
  %
  % returns a cell of selections
  
  s = {}; this.event = '';
  while ~isempty(this) && isstruct(this)

    if strcmpi(this.event, 'rotate')
      r=rotate3d(gca);
      disp([ mfilename ': switching to rotate mode. Deselect that mode (rotate icon or Tools menu) to continue selecting the ROI.' ]);
      set(r, 'Enable','on');
      waitfor(r, 'Enable','off');
    end
    if ~isempty(this)
      disp([ mfilename ': adding new polygon ' num2str(numel(s)+1) ]);
      s{end+1} = this; 
    end
    if strcmpi(this.event, 'end')
      return
    end
    this = selectpoints_ginput();
  end % while
end % selectpoints
  
function s = selectpoints_ginput()
% selectpoints: define points on the current axis.
%   mouse left-click                add a point/line
%   right-click or DEL or BACKSPACE removes last point
%   return or 'q'                   terminate input
%   ESC                             abort
%
% returns:
%   a structure with members XV,YV,ZV containing the vertices
  button   = 1;
  s.XV = []; s.YV = []; s.ZV = []; s.h = []; s.event = '';

  while button
    % handle KeyPressFcn events
    if ~isfield(s, 'XV'),   return; end  % quit when e.g. ESC pressed
    
    % we collect points on the axis and draw lines
    [x, y, button] = ginput(1);
    
    % handle the returned data
    % button can be 1-3 for the mouse or a pressed key as integer
    if isempty(x) || isempty(button)  % hit Return
      s.event = 'end';
      return
    end
    switch lower(button)
    case { 3, sprintf('\b') }       % BS  back/erase
      if ~isempty(s.XV)
        s.XV(:,end) = [];
        s.YV(:,end) = [];
        if ~isempty(s.ZV) s.ZV(:,end) = []; end
        delete(s.h(end)); s.h(end) = [];
      end
      continue
    case {'c', 127 }                % DEL or 'c' erase all
      delete(s.h);
      s.XV = []; s.YV = []; s.ZV = []; s.h = [];
      continue
    case { 2, 'q', 'x', sprintf('\r') }     % CR  return or 'q' end
      s.event = 'end';
      return
    case 27                         % ESC escape/abort
      % quit -> s = [];
      delete(s.h);
      s = []; return
    case 'r'
      % toggle rotate tool
      s.event = 'rotate';
      return
    case 'a'
      % go to a next polygon
      s.event = 'next';
      return
    case 'h'                        % help
      helpdlg({[ mfilename ': Add vertices (points) on the axis to', ...
      'define a Region Of Interest area' ], ...
      ' ', ...
      '  mouse left-click           add a point/line', ...
      '  right-click or BACKSPACE   removes last point', ...
      '  DEL or "c"                 removes all points (clear)', ...
      '  return or "q" or middle-click: terminate input (quit)', ...
      '  a                          define a new separate polygon (add)', ...
      '  r                          change orientation (2D/3D). Disable rotation to continue ROI.', ...
      '  ESC                        abort'
      }, [ mfilename ': Define ROI' ]);
      continue
    end

    % store the point
    cc = get(gca,'CurrentPoint'); % 2 rows
    xv  = cc(:,1); yv=cc(:,2);
    % draw a line
    if size(cc, 2) == 2 % 2D axis (flat)
      h = line(xv,yv); zv = [];
    else                % 3D axis
      zv = cc(:,3);
      h = line(xv,yv,zv);
    end
    set(h, 'Color','r', 'Marker', 's', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
    
    % we add these vertices so that later we search data points in the polygon
    s.XV = [ s.XV xv ];
    s.YV = [ s.YV yv ];
    s.ZV = [ s.ZV zv ];
    s.h  = [ s.h  h  ];
    s.event = 'point';
  end
end % selectpoints_ginput

% ------------------------------------------------------------------------------
function hc = get_children(h)
% get_children: find all children objects that have xdata/ydata properties and return their ID
  
  hc = [];
  
  for index=1:numel(h)
    try
      get(h(index), 'XData');
      hc = [ hc h(index) ];
    catch
      % skip that one
    end
    % look for any Children
    hc = [ hc get_children(get(h(index),'Children')) ];
  end
end % get_children
 
% ------------------------------------------------------------------------------
% https://fr.mathworks.com/matlabcentral/fileexchange/10226-inhull
% BSD license
function in = inhull(testpts,xyz,tess,tol)

% inhull: tests if a set of points are inside a convex hull
% usage: in = inhull(testpts,xyz)
% usage: in = inhull(testpts,xyz,tess)
% usage: in = inhull(testpts,xyz,tess,tol)
%
% arguments: (input)
%  testpts - nxp array to test, n data points, in p dimensions
%       If you have many points to test, it is most efficient to
%       call this function once with the entire set.
%
%  xyz - mxp array of vertices of the convex hull, as used by
%       convhulln.
%
%  tess - tessellation (or triangulation) generated by convhulln
%       If tess is left empty or not supplied, then it will be
%       generated.
%
%  tol - (OPTIONAL) tolerance on the tests for inclusion in the
%       convex hull. You can think of tol as the distance a point
%       may possibly lie outside the hull, and still be perceived
%       as on the surface of the hull. Because of numerical slop
%       nothing can ever be done exactly here. I might guess a
%       semi-intelligent value of tol to be
%
%         tol = 1.e-13*mean(abs(xyz(:)))
%
%       In higher dimensions, the numerical issues of floating
%       point arithmetic will probably suggest a larger value
%       of tol.
%
%       DEFAULT: tol = 0
%
% arguments: (output)
%  in  - nx1 logical vector
%        in(i) == 1 --> the i'th point was inside the convex hull.
%  
% Example usage: The first point should be inside, the second out
%
%  xy = randn(20,2);
%  tess = convhulln(xy);
%  testpoints = [ 0 0; 10 10];
%  in = inhull(testpoints,xy,tess)
%
% in = 
%      1
%      0
%
% A non-zero count of the number of degenerate simplexes in the hull
% will generate a warning (in 4 or more dimensions.) This warning
% may be disabled off with the command:
%
%   warning('off','inhull:degeneracy')
%
% See also: convhull, convhulln, delaunay, delaunayn, tsearch, tsearchn
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 3.0
% Release date: 10/26/06

% get array sizes
% m points, p dimensions
  p = size(xyz,2);
  [n,c] = size(testpts);
  if p ~= c
    error 'testpts and xyz must have the same number of columns'
  end
  if p < 2
    error 'Points must lie in at least a 2-d space.'
  end

  % was the convex hull supplied?
  if (nargin<3) || isempty(tess)
    tess = convhulln(xyz);
  end
  [nt,c] = size(tess);
  if c ~= p
    error 'tess array is incompatible with a dimension p space'
  end

  % was tol supplied?
  if (nargin<4) || isempty(tol)
    tol = 0;
  end

  % build normal vectors
  switch p
    case 2
      % really simple for 2-d
      nrmls = (xyz(tess(:,1),:) - xyz(tess(:,2),:)) * [0 1;-1 0];
      
      % Any degenerate edges?
      del = sqrt(sum(nrmls.^2,2));
      degenflag = (del<(max(del)*10*eps));
      if sum(degenflag)>0
        warning('inhull:degeneracy',[num2str(sum(degenflag)), ...
          ' degenerate edges identified in the convex hull'])
        
        % we need to delete those degenerate normal vectors
        nrmls(degenflag,:) = [];
        nt = size(nrmls,1);
      end
    case 3
      % use vectorized cross product for 3-d
      ab = xyz(tess(:,1),:) - xyz(tess(:,2),:);
      ac = xyz(tess(:,1),:) - xyz(tess(:,3),:);
      nrmls = cross(ab,ac,2);
      degenflag = false(nt,1);
    otherwise
      % slightly more work in higher dimensions, 
      nrmls = zeros(nt,p);
      degenflag = false(nt,1);
      for i = 1:nt
        % just in case of a degeneracy
        % Note that bsxfun COULD be used in this line, but I have chosen to
        % not do so to maintain compatibility. This code is still used by
        % users of older releases.
        %  nullsp = null(bsxfun(@minus,xyz(tess(i,2:end),:),xyz(tess(i,1),:)))';
        nullsp = null(xyz(tess(i,2:end),:) - repmat(xyz(tess(i,1),:),p-1,1))';
        if size(nullsp,1)>1
          degenflag(i) = true;
          nrmls(i,:) = NaN;
        else
          nrmls(i,:) = nullsp;
        end
      end
      if sum(degenflag)>0
        warning('inhull:degeneracy',[num2str(sum(degenflag)), ...
          ' degenerate simplexes identified in the convex hull'])
        
        % we need to delete those degenerate normal vectors
        nrmls(degenflag,:) = [];
        nt = size(nrmls,1);
      end
  end

  % scale normal vectors to unit length
  nrmllen = sqrt(sum(nrmls.^2,2));
  % again, bsxfun COULD be employed here...
  %  nrmls = bsxfun(@times,nrmls,1./nrmllen);
  nrmls = nrmls.*repmat(1./nrmllen,1,p);

  % center point in the hull
  center = mean(xyz,1);

  % any point in the plane of each simplex in the convex hull
  a = xyz(tess(~degenflag,1),:);

  % ensure the normals are pointing inwards
  % this line too could employ bsxfun...
  %  dp = sum(bsxfun(@minus,center,a).*nrmls,2);
  dp = sum((repmat(center,nt,1) - a).*nrmls,2);
  k = dp<0;
  nrmls(k,:) = -nrmls(k,:);

  % We want to test if:  dot((x - a),N) >= 0
  % If so for all faces of the hull, then x is inside
  % the hull. Change this to dot(x,N) >= dot(a,N)
  aN = sum(nrmls.*a,2);

  % test, be careful in case there are many points
  in = false(n,1);

  % if n is too large, we need to worry about the
  % dot product grabbing huge chunks of memory.
  memblock = 1e6;
  blocks = max(1,floor(n/(memblock/nt)));
  aNr = repmat(aN,1,length(1:blocks:n));
  for i = 1:blocks
     j = i:blocks:n;
     if size(aNr,2) ~= length(j),
        aNr = repmat(aN,1,length(j));
     end
     in(j) = all((nrmls*testpts(j,:)' - aNr) >= -tol,1)';
  end

end % inhull
