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
      
      % Matlab <= 2014a:
      % errorbar returns a hggroup. Its children are 2 'Line' objects.
      
      % Matlab >= 2014b:
      % with new HandleGraphics (HG2), errorbars are single objects, not two separate
      % lines. To be able to hide the errorbar, we also plot the line alone.
      
      % we first call errorbar. If the result is a not an hggroup, but a single 
      % Errorbar object, we set its LineStyle to None and plot separately the line.
      % then we create an hggroup
      h=[];
      try
        if ~isempty(this_method)
          h = errorbar(x,y,e,this_method, varargin{:});
        end
      end
      if isempty(h) % if errorbar(this_method) failed or no 'this_method'
        h = errorbar(x,y,e, varargin{:});
      end
      % handle the HG2 case
      if ishandle(h) && strcmp(get(h,'Type'), 'errorbar')
        % HG2 errorbar
        set(h,'LineStyle','none');  % we only keep the error bars
        % we plot the line separately
        h2=[];
        hold on
        try
          if ~isempty(this_method)
            h2 = plot(x,y,this_method, varargin{:});
          end
        end
        if isempty(h2) % plot(this_method) failed or no 'this_method'
          h2 = plot(x,y, varargin{:});
        end
        if ishandle(h2)
          h = [h2 h];
          hg = hggroup;
          set(h, 'Parent', hg);
          h = hg;
        end
      end % HG2 case
      
      % option to hide errorbar
      if ~isempty(strfind(method, 'hide_err')) || all(abs(e) >= abs(y) | e == 0) % && numel(h) > 1
        eh = findobj(h,'Type','errorbar');
        if isempty(eh) && numel(h) > 1, eh = h(2);
        elseif isempty(eh) && strcmp(get(h,'Type'),'hggroup')
          eh = get(h, 'Children');
          if numel(eh) > 1, eh = eh(2); end
        end
        set(eh, 'Visible','off');
      end
    end
  end
  
  % add manual ticks when specified
  if isfield(a, 'XTickLabel')
    xtick = get(a,'XTickLabel');
    if ischar(xtick) && isvector(xtick)
      xtick = textscan(xtick, '%s','delimiter', ' ');
      xtick = xtick{1};
    end
    set(gca, 'XTickLabel', xtick);
  end
