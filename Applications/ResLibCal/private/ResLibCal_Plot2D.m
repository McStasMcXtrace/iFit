function out=ResLibCal_Plot2D(out, mode)
%
% MATLAB function to plot the projections of the resolution ellipse
% of a triple axis
%
% Input:
%  out:  EXP ResLib structure 
%  mode: can be set to 'rlu' so that the plot is in lattice RLU frame [abc] instead of [xyz]
%        can also specify 'qz' to indicate x,y,z in place of x,y,E mode.

% input parameters
if nargin < 1, out = []; end
if nargin < 2, mode=''; end
if isempty(out), return; end

if ~isfield(out, 'resolution') 
  disp([ mfilename ': Input argument does not contain resolution. Skipping.' ]); 
  return
end
if isfield(out, 'EXP');
  EXP = out.EXP;
else
  out
  disp([ mfilename ': Input argument does not contain EXP structure. Skipping.' ])
  return
end

% handle resolution for scans
if numel(out.resolution) == 1
  resolutions = { out.resolution };
else
  resolutions = out.resolution;
end

max_points = ceil(max(100, 300/numel(resolutions))); % nb of MC points per cloud in 3 panes

for index=1:numel(resolutions)
  resolution = resolutions{index};
  
  if ~resolution.R0, continue; end
  H=resolution.HKLE(1); K=resolution.HKLE(2); 
  L=resolution.HKLE(3); W=resolution.HKLE(4);

  if ~isempty(strfind(mode,'rlu'))
    frame = resolution.rlu;
  else
    frame = resolution.spec;
  end
  NP       = frame.RM;
  FrameStr = {frame.frameStr{:},'\omega meV'};
  [Labels, FrameStr] = strtok(FrameStr);
  cloud    = frame.cloud;
  Units    = frame.unit;
  centre   = frame.Q; centre(4) = W;
  FrameStr{end+1} = ''; % no specific axis for energy (4)
  if index > 1, Labels = []; end

  if isempty(NP) || ~all(isreal(NP)), return; end
  
  if isempty(strfind(mode,'cloud')), cloud=[]; end
  
  % plot the 3 subplots for projections
  % each plot is shown as a gauss2d, contour.
  % add context menus accordingly
  % arguments: Subpanel index, NP(ix,iy), 
  [dummy, NP2] = rc_int(4,1, NP); % this function strips out the row=col=4, and corrects determinant
                                  % using the cofactor rule.
  [dummy, NP2] = rc_int(3,1, NP2);
  ResLibCal_Proj_plot2D(1, 1,2, NP2, Labels, FrameStr, Units, 'xy', cloud, centre, max_points);  % xy
  if index == 1, title([ 'ResLibCal ' datestr(now) ]); end
  if index < numel(resolutions), hold on; else hold off; end
  % if numel(resolutions) > 1 && index==1, colorbar('North'); end
  
  if ~isempty(strfind(mode,'qz'))
    [dummy, NP2] = rc_int(4,1, NP);
    [dummy, NP2] = rc_int(2,1, NP2);
    ResLibCal_Proj_plot2D(2, 1,3, NP2, Labels, FrameStr, Units, 'xz', cloud, centre, max_points);  % xe
    if index == 1, title('Vertical Q_z resolution'); end
    if index < numel(resolutions), hold on; else hold off; end
    
    [dummy, NP2] = rc_int(4,1, NP);
    [dummy, NP2] = rc_int(1,1, NP2);
    ResLibCal_Proj_plot2D(3, 2,3, NP2, Labels, FrameStr, Units, 'yz', cloud, centre, max_points);  % ye
  else
    [dummy, NP2] = rc_int(3,1, NP);
    [dummy, NP2] = rc_int(2,1, NP2);
    ResLibCal_Proj_plot2D(2, 1,4, NP2, Labels, FrameStr, Units, 'xz', cloud, centre, max_points);  % xe
    if index == 1, title('Energy resolution'); end
    if index < numel(resolutions), hold on; else hold off; end
    
    [dummy, NP2] = rc_int(3,1, NP);
    [dummy, NP2] = rc_int(1,1, NP2);
    ResLibCal_Proj_plot2D(3, 2,4, NP2, Labels, FrameStr, Units, 'yz', cloud, centre, max_points);  % ye
  end
  if index == 1, title(EXP.method); end
  if index < numel(resolutions), hold on; else hold off; end

end % for
hold off;

% now add the text box
% display a text edit uicontrol so that users can select/copy/paste
[res, inst] = ResLibCal_FormatString(out, mode);
message = [ res; inst ];

% fill 4th sub-panel with uicontrol
p(1) = 0.5; p(2) = 0.01; p(3) = 0.49; p(4) = 0.49;
h = findobj(gcf, 'Tag','ResLibCal_View2_Edit');
if isempty(h)
  h = uicontrol('Tag','ResLibCal_View2_Edit', ...
    'Style','edit','Units','normalized',...
    'Position',p, 'Max',2,'Min',0, ...
    'String', message, 'FontName', 'FixedWidth', ...
    'FontSize',8,'BackgroundColor','white','HorizontalAlignment','left');
  % contextual menu
  uicm = uicontextmenu;
  uimenu(uicm, 'Label', 'ResLibCal: Resolution: results') ;
  uimenu(uicm, 'Label', 'Copy to clipboard', ...
    'Callback', 'tmp_s=get(gco,''String''); clipboard(''copy'', sprintf(''%s\n'', tmp_s{:})); clear tmp_s;');
  set(h, 'UIContextMenu', uicm, 'Tag','ResLibCal_Proj_Context_text');
else
  set(h, 'String',message);
end

% ------------------------------------------------------------------------------
function h=ResLibCal_Proj_plot2D(isub, ix,iy,NP, Labels, FrameStr, Units, panel_name, cloud, centre, max_points)
% plot the subpanel isub, with resolution projection(ix,iy)
% each plot is a gaussian 2D, with axes FrameStr{ix|iy}
% expressed in Units, and short indices name as 'panel_name'

  persistent ResLibCal_version
  
  if isempty(ResLibCal_version), ResLibCal_version = ResLibCal('version'); end
  
  % reduce resolution matrix to (ix,iy)
  % and plot
  subplot(2,2,isub);
  if ~isempty(Labels); cla; end
  % plot the cloud projection if any
  if ~isempty(cloud)
    x=cloud{ix}; y=cloud{iy};
    if isempty(find([ix iy] == 4)), e = cloud{4};
    else e = cloud{3}; end
    if numel(x) > max_points, x=x(1:max_points); end
    if numel(y) > max_points, y=y(1:max_points); end
    if numel(e) > max_points, e=e(1:max_points); end
    h=scatter(x,y,3, e);
    set(h,'DisplayName','cloud');
    hold on
  end
  % plot the ellipse on top.
  [h, XX] = Ellipse_plot(NP, centre([ix iy]));
  set(h,'DisplayName',panel_name);
  if ~isempty(Labels)
  
    X=XX(1,:); sx=max(X(:)) - min(X(:));
    Y=XX(2,:); sy=max(Y(:)) - min(Y(:));

    xlabel({ [ '{\bf ' Labels{ix} '} ' FrameStr{ix} ' [' Units ']' ], [ '{\delta}' Labels{ix} '=' num2str(sx, 3) ]})
    ylabel({ [ '{\bf ' Labels{iy} '} ' FrameStr{iy} ' [' Units ']' ], [ '{\delta}' Labels{iy} '=' num2str(sy, 3) ]})

    grid on

    % contextual menu
    if isempty(findobj(gca, 'Tag',[ 'ResLibCal_Proj_Context_' panel_name ]))
      uicm = uicontextmenu;
      uimenu(uicm, 'Label', [ 'ResLibCal: Resolution: Q' panel_name(1) 'Q' panel_name(2) ' projection' ]) ;
      uimenu(uicm, 'Separator','on', 'Label', 'Duplicate View...', 'Callback', ...
         [ 'tmp_cb.g=gca;' ...
           'tmp_cb.f=figure; tmp_cb.c=copyobj(tmp_cb.g,gcf); ' ...
           'set(tmp_cb.c,''position'',[ 0.1 0.1 0.85 0.8]);' ...
           'set(gcf,''Name'',''Copy of ResLibCal: ' panel_name ' Resolution ' Units '''); ' ...
           'set(gca,''XTickLabelMode'',''auto'',''XTickMode'',''auto'');' ...
           'set(gca,''YTickLabelMode'',''auto'',''YTickMode'',''auto'');' ...
           'set(gca,''ZTickLabelMode'',''auto'',''ZTickMode'',''auto'');']);
      uimenu(uicm, 'Label','Toggle grid', 'Callback','grid');
      uimenu(uicm, 'Label','Reset View', 'Callback','view(2);alpha(0.5);axis tight;rotate3d off;');
      uimenu(uicm, 'Separator','on','Label', 'About ResLibCal...', ...
        'Callback',[ 'msgbox(''' ResLibCal_version ''',''About ResLibCal'',''help'')' ]);
      set(gca, 'UIContextMenu', uicm, 'Tag',[ 'ResLibCal_Proj_Context_' panel_name ]);
    end
  end
  
  
% ------------------------------------------------------------------------------
 
