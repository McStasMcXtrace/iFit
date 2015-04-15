function [h, xlab, ylab, ret] = iData_plot_1d(a, method, this_method, varargin)
% iData_plot_1d: plot a 1D iData object
% used in iData/plot

  ret = 0;

  if size(a,1) ==1 && size(a,2) > 1
    a = transpose(a);
  end
  [x, xlab] = getaxis(a,1); x=double(x(:));
  [y, ylab] = getaxis(a,0); y=double(y(:));
  e         = get(a,'Error');   e=real(double(e(:)));
  m         = get(a,'Monitor'); m=real(double(m(:)));
  if not(all(m == 1 | m == 0)),
    e=genop(@rdivide,e,m); ylab = [ylab ' per monitor' ];
  end
  y=real(y);
  y(isinf(y)) = nan;
  
  if isempty(method), method='b-'; end
  % handle side-by-side 1D plots
  if ~isempty(strfind(method,'plot3'))    || ~isempty(strfind(method,'stem')) ...
   || ~isempty(strfind(method,'scatter')) || ~isempty(strfind(method,'mesh')) ...
   || ~isempty(strfind(method,'surf') )   || ~isempty(strfind(method,'waterfall'))
  	ax = getaxis(a,2);
  	if isempty(ax)
  		ax = 0;
    end
    if length(ax) == 1
    	ax = ax*ones(size(a));
    end
    % need to create this axis
    setalias(a, 'Axis_2', ax);
    setaxis(a, 2, 'Axis_2');
    h = plot(a, method, varargin{:});
    ret = 1;
  else
    if all(e == 0 | ~isfinite(e)) || length(x) ~= length(e)
      % no errorbar -> single line
      if length(this_method)
        try
          h = plot(x,y, this_method, varargin{:});
        catch
          this_method=[];
        end
      end
      if ~length(this_method) h = plot(x,y, varargin{:}); end
    else
      % with errorbar -> plot both with and without errorbar so that we can
      % hide the errorbar, retaining the data set only

      % the returned handle is the plot, and then any errorbar stuff (hidable)
      hg = hggroup; h=[];
      if length(this_method)
        try
          h = [ plot(x,y,this_method, varargin{:}, 'Parent',hg) ...
                errorbar(x,y,e,this_method, varargin{:}, 'Parent',hg) ];
        catch
          delete(h)
          this_method=[]; % indicate we failed so that we try without line option
        end
      end
      if ~length(this_method)
        try
          h = [ plot(x,y,varargin{:}, 'Parent',hg) ...
                errorbar(x,y,e,varargin{:}, 'Parent',hg) ];
        catch ME
          delete(h)
          h = errorbar(x,y,e, varargin{:});
        end
      end
      
      % make sure the linespec of the single plot and errorbar coincide
      if numel(h) > 1 
        l = findobj(h,'Type','line');
        er = findobj(h,'Type','errorbar');
        specs = {'LineStyle','','LineWidth','','Color', '', ...
          'MarkerEdgeColor','','MarkerFaceColor','','MarkerSize','' };
        specs_values = specs;
        for index=1:numel(specs)
          if ~isempty(specs{index})
            specs_values{index+1} = get(l(1), specs{index});
          end
        end
        set(er, specs_values{:});
      end
      
      % option to hide errorbar
      if (~isempty(strfind(method, 'hide_err')) || all(abs(e) >= abs(y) | e == 0)) && numel(h) > 1
        set(findobj(h,'Type','errorbar'), 'Visible','off');
      end
    end
  end
