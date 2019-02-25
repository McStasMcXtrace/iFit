function s = read_fig(filename)
% read_fig Wrapper to directly read Matlab Figures
%  s = read_fig(filename)
%
% Input:  filename: FIG filename or graphic handle
% output: structure
% Example: h=surf(peaks); s=read_fig(h); close(gcf); isstruct(s)
%
% (c) E.Farhi, ILL. License: EUPL.
% See also: read_idl, read_lvm, read_tdms, read_igor

s=[];

if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
    Matlab_FIG.name     ='Matlab Figure';
    Matlab_FIG.method   =mfilename;
    Matlab_FIG.extension='fig';
    Matlab_FIG.postprocess='load_fig';
    s=Matlab_FIG;
    return
end

if ishandle(filename)
  s = handle2iData(filename);
else
  f       = openfig(filename, 'new','invisible');
  s = handle2iData(f);
end

% ------------------------------------------------------------------------------
% handle2iData: converts a figure/graphic handle into a structure
function out=handle2iData(in)
  try 
    t = get(in,'DisplayName');
    if isempty(t), t=get(get(in,'Parent'),'DisplayName'); end
  catch
    t=[]; end
  if isempty(t), t=get(in,'Tag'); end
  if isempty(t), 
      t=num2str(double(in)); 
  end
  if strcmp(get(in,'type'),'hggroup')
    t = [ 'figure ' t ];
    h = get(in,'Children');
    out = handle2iData(h(1)); % first item
    out.Title=t;
    out.Label=t;
  elseif strcmp(get(in,'type'),'line')
    x = get(in,'xdata'); 
    y = get(in,'ydata'); 
    index = find(~isnan(x) & ~isnan(y));
    if length(index)~=numel(x), x = x(index); y=y(index); end
    c = get(in,'color');
    m = get(in,'marker');
    l = get(in,'linestyle');
    out.x=x; out.Signal=y;
    try xl = get(get(in,'parent'),'XLabel'); xl=get(xl,'String'); catch 
        xl='x'; end; 
    try yl = get(get(in,'parent'),'YLabel'); yl=[ get(yl,'String') ' ' ]; catch 
        yl=''; end;
    try tl = get(get(in,'parent'),'Title');  tl=[ get(tl,'String') ' ' ]; catch 
        tl=''; end;
    out.xlabel=xl; out.ylabel=yl;
    t = [ 'line ' t ];
    out.Title = [ tl yl t ];
    out.DisplayName = t;
    out.Label=[ t ' marker ' m ' color ' num2str(c) ];
  elseif strcmp(get(in,'type'),'image')
    x = get(in,'xdata'); 
    y = get(in,'ydata');
    z = get(in,'cdata');
    t = [ 'image ' t ];
    out.x=x; out.y=y; out.Signal=z;
    try xl = get(get(in,'parent'),'XLabel'); xl=get(xl,'String'); catch 
        xl='x'; end
    try yl = get(get(in,'parent'),'YLabel'); yl=get(yl,'String'); catch 
        yl='y'; end
    try zl = get(get(in,'parent'),'ZLabel'); zl=[ get(zl,'String') ' ' ]; catch 
        zl=''; end 
    try tl = get(get(in,'parent'),'Title');  tl=[ get(tl,'String') ' ' ]; catch 
        tl=''; end
    out.xlabel=xl; out.ylabel=yl; out.zlabel=zl;
    out.Title = t;
    out.DisplayName = tl;
    out.Label=t;
  elseif strcmp(get(in,'type'),'surface')
    x = get(in,'xdata'); 
    y = get(in,'ydata'); 
    z = get(in,'zdata'); 
    c = get(in,'cdata'); 
    % index=find(~isnan(x) & ~isnan(y) & ~isnan(z) & ~isnan(c)); 
    % if length(index)~=prod(size(x)), x = x(index); y=y(index); z=z(index); c=c(index); end
    l = get(in,'linestyle');
    if all(z == c)
      out.x=x; out.y=y; out.Signal=z;
    else
      out.x=x; out.y=y; out.z=z; out.Signal=c;
    end
    try 
        xl = get(get(in,'parent'),'XLabel'); xl=get(xl,'String'); 
    catch 
        xl='x'; end
    try 
        yl = get(get(in,'parent'),'YLabel'); yl=get(yl,'String'); 
    catch
        yl='y'; end
    try 
        zl = get(get(in,'parent'),'ZLabel'); zl=[ get(zl,'String') ' ' ]; 
    catch 
        zl=''; end 
    try 
        tl = get(get(in,'parent'),'Title');  tl=[ get(tl,'String') ' ' ]; 
    catch 
        tl=''; end

    out.xlabel=xl; out.ylabel=yl; out.zlabel=zl;
    if all(z == c)
      t = [ tl zl t ];
    else
      if isempty(zl), zl='z'; end
      zlabel(out, zl);
      t = [ tl t ];
    end
    t = [ 'surface ' t ];
    out.Title = t;
    out.DisplayName = tl;
    out.Label=[ t ' line ' l ];
  else
    h = [ findobj(in, 'type','line') ; findobj(in, 'type','surface') ; findobj(in, 'type','image')  ];
    out = [];
    for index=1:length(h)
      this_out = handle2iData(h(index));
      if isempty(this_out.Title) && ~isempty(t)
        this_out.Title = t;
        this_out.Label = t;
        this_out.DisplayName = t;
      end
      if ~isempty(t), this_out.Source = t; end
      if  ~isscalar(get(this_out,'Signal'))
        out = [ out this_out ];
      end
    end
  end
% ------------------------------------------------------------------------------
