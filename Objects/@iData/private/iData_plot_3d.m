function [h, xlab, ylab, zlab, ret] = iData_plot_3d(a, method, this_method, varargin)
% iData_plot_3d: plot a 3D iData object
% used in iData/plot

ret = 0; 
if ndims(a) == 1
  zlab='';
  if strcmp(method, 'scatter'), method='plot'; end
  [h, xlab, ylab, ret] = iData_plot_1d(a, method, this_method, varargin{:}); % in private
  return
elseif ndims(a) == 2
  [h, xlab, ylab, zlab] = iData_plot_2d(a, method, this_method, varargin{:}); % in private
  return
end

% first test if this is an image
  if isfield(a.Data,'cdata')
    h=image(a.Data.cdata);
    xlab=''; ylab=''; zlab='';
  else
    % check if a rebining on a grid is required
    if (~isvector(a) && (~isempty(strfind(method,'plot3')) || ~isempty(strfind(method,'scatter')) ))
      a = meshgrid(a); % make sure we get a grid
    end
    [x, xlab] = getaxis(a,2); x=double(x);
    [y, ylab] = getaxis(a,1); y=double(y);
    [z, zlab] = getaxis(a,3); z=double(z);
    [c, clab] = getaxis(a,0); c=double(c); % c(isinf(c)) = nan;
    m         = get(a,'Monitor');
    if not(all(m(:) == 1 | m(:) == 0)), clab = [clab ' per monitor' ]; end
    if isvector(a) >= 3 || ~isempty(strfind(method, 'scatter')) % plot3-like
      if isempty(strfind(method, 'plot3'))
        h = hggroup;
        h3=fscatter3(x(:),y(:),z(:),c(:), this_method);     % scatter3: may require meshgrid
        set(h3,'Parent',h);
      else
        h=plot3(x(:),y(:),z(:), this_method, varargin{:});
      end
    else
      if ~isempty(strfind(method, 'plot3')) % vol3d: does not require meshgrid
        h = hggroup;
        h3 = vol3d('cdata',c,'texture','3D','xdata',x,'ydata',y,'zdata',z);
        alphamap('vdown'); % make object transparent on borders and solid in center
        h3 = vol3d(h3);
        set(h3.handles,'Parent',h);
      elseif ~isempty(strfind(method, 'waterfall')) || ~isempty(strfind(method, 'contour'))
        h = hggroup;
        if ~isempty(strfind(method, ' y '))
          iy = linspace(min(y(:)), max(y(:)), 10);
          hc = contourslice(x,y,z,c,[],iy,[], varargin{:});
        elseif ~isempty(strfind(method, ' x '))
          ix = linspace(min(x(:)), max(x(:)), 10);
          hc = contourslice(x,y,z,c,ix,[],[], varargin{:});
        else
          iz = linspace(min(z(:)), max(z(:)), 10);
          hc = contourslice(x,y,z,c,[],[],iz, varargin{:});
        end
        set(hc,'Parent',h);
      elseif ~isempty(strfind(method, 'slice')) % sliceomatic
        slice(a); h=[];
      else % method='surf'
        % isosurface: require meshgrid
        if ~isempty(strfind(method, 'mean'))
          iso = mean(c(:));
        elseif ~isempty(strfind(method, 'half'))
          iso = (min(c(:))+max(c(:)))/2;
        elseif ~isempty(strfind(method, 'median'))
          iso = median(c(:));
        else
          iso = [];
        end
        try
          if ~isempty(iso), 
            isosurface(x,y,z, c, iso, varargin{:});
          else 
            isosurface(x,y,z, c, varargin{:}); 
          end
        end
        h = findobj(gca,'type','patch');
      end
    end
    view(3);
    zlabel(zlab);
  end
