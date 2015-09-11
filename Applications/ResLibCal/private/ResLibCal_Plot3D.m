function out = ResLibCal_Plot3D(out, mode)
% adapted from ResLib
%===================================================================================
%  function ResPlot3D(
%         H,K,L,W,EXP,RANGE,EllipsoidStyle,XYStyle,XEStyle,YEStyle,SMA,SMAp,SXg,SYg)
%  ResLib v.3.4
%===================================================================================
%
% For a specified scan, plot pretty  resolution ellipsoids in 3D.
% If a SMA cross section is specified, the calculated dispersion is plotted as well.
%
% A. Zheludev, 1999-2007
% Oak Ridge National Laboratory
%====================================================================================

% Calls: StandardSystem, GetLattice, scalar, modvec
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
  disp([ mfilename ': Input argument does not contain EXP structure. Skipping.' ])
  return
end

% clean up current axis if redraw
if isempty(strfind(mode,'scan')) && ~isempty(findobj(gcf,'Tag','ResLibCal_View3_Context'))
  delete(findobj(gcf,'Tag','ResLibCal_View3_Volume'));
  delete(findobj(gcf,'Tag','ResLibCal_View3_Cloud'));
end

% handle resolution for scans
if numel(out.resolution) == 1
  resolutions = { out.resolution };
else
  resolutions = out.resolution;
end

cla;

max_points = max(100, 1000/numel(resolutions)); % nb of MC points per cloud

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
  
  % reduce dimensionality of the resolution mtrix
  if ~isempty(strfind(mode,'qz'))
    [dummy, NP] = rc_int(4,1, NP); % this function strips out the row=col=4, and corrects determinant
                                   % using the cofactor rule.
    ResLibCal_Proj_plot3D([1 2 3], NP, Labels, FrameStr, Units, cloud, centre, max_points);  % xyz
    if index == 1, title([ 'Resolution in ' Labels{1:3} ' - ' out.EXP.method ]); end
    if index < numel(resolutions), hold on; else hold off; end
  else
    [dummy, NP] = rc_int(3,1, NP);
    ResLibCal_Proj_plot3D([1 2 4], NP, Labels, FrameStr, Units, cloud, centre, max_points);  % xye
    if index == 1, title([ 'Resolution in ' Labels{[1 2 4]} ' - ' out.EXP.method ]); end
    if index < numel(resolutions), hold on; else hold off; end
  end
  
end % for

% ------------------------------------------------------------------------------
function h=ResLibCal_Proj_plot3D(index, NP, Labels, FrameStr, Units, cloud, centre, max_points)

  if ~isempty(Labels), cla; end
  % plot the cloud projection if any
  if ~isempty(cloud)
    x=cloud{index(1)}; y=cloud{index(2)}; z=cloud{index(3)};
    if isempty(find(index == 4)), e = cloud{4};
    else e = cloud{3}; end
    if numel(x) > max_points, x=x(1:max_points); end
    if numel(y) > max_points, y=y(1:max_points); end
    if numel(z) > max_points, z=z(1:max_points); end
    if numel(e) > max_points, e=e(1:max_points); end
    h=scatter3(x,y,z,3,e,'MarkerFaceColor','b');
    % plot3(x,y,z,'o');
    if ~isempty(Labels)
      set(h,'DisplayName',[ Labels{index} ' (cloud)' ], 'Tag', 'ResLibCal_View3_Cloud');
    end
    hold on
  end
  
  % plot the ellipse on top.
  [h, XX,YY,ZZ] = Ellipse_plot(NP, centre(index),10);
  if ~isempty(Labels)
    set(h,'DisplayName',[ Labels{index} ],'Tag','ResLibCal_View3_Volume');
    
    sx=max(XX(:)) - min(XX(:)); ix=index(1);
    sy=max(YY(:)) - min(YY(:)); iy=index(2);
    sz=max(ZZ(:)) - min(ZZ(:)); iz=index(3);
    
    xlabel({[ '{\bf ' Labels{ix} '} ' FrameStr{ix} ' [' Units ']' ], ...
         [ '{\delta}' Labels{ix} '=' num2str(sx, 3) ]})
    ylabel({[ '{\bf ' Labels{iy} '} ' FrameStr{iy} ' [' Units ']' ], ...
         [ '{\delta}' Labels{iy} '=' num2str(sy, 3) ]})
    zlabel({[ '{\bf ' Labels{iz} '} ' FrameStr{iz} ' [' Units ']' ], ...
         [ '{\delta}' Labels{iz} '=' num2str(sz, 3) ]})
         
    grid on
    
    % add contextual menu
    if isempty(findobj(gcf,'Tag','ResLibCal_View3_Context'))
      %finalize 3D plot
      box on; grid on;
      view(3);

      uicm = uicontextmenu;
      uimenu(uicm, 'Label', 'ResLibCal: Resolution: 3D') ;
      uimenu(uicm, 'Separator','on', 'Label', 'Duplicate View...', 'Callback', ...
         [ 'tmp_cb.g=gca;' ...
           'tmp_cb.f=figure; tmp_cb.c=copyobj(tmp_cb.g,gcf); ' ...
           'set(tmp_cb.c,''position'',[ 0.1 0.1 0.85 0.8]);' ...
           'set(gcf,''Name'',''Copy of ResLibCal: TAS view''); ' ...
           'set(gca,''XTickLabelMode'',''auto'',''XTickMode'',''auto'');' ...
           'set(gca,''YTickLabelMode'',''auto'',''YTickMode'',''auto'');' ...
           'set(gca,''ZTickLabelMode'',''auto'',''ZTickMode'',''auto'');']);
      uimenu(uicm, 'Label','Toggle grid', 'Callback','grid');
      uimenu(uicm, 'Label','Reset Flat/3D View', 'Callback', [ ...
          '[tmp_a,tmp_e]=view; if (tmp_a==0 & tmp_e==90) view(3); else view(2); end;' ...
          'clear tmp_a tmp_e; lighting none; shading flat;' ]);
      uimenu(uicm, 'Label','Smooth View','Callback', 'shading interp;');
      uimenu(uicm, 'Label','Add Light','Callback', 'light;lighting phong;');
      uimenu(uicm, 'Label','Add Transparency','Callback', 'alphamap(''decrease''); for tmp_h=get(gca, ''children'')''; try; alpha(tmp_h,0.7*get(tmp_h, ''facealpha'')); end; end;');
      uimenu(uicm, 'Label','Toggle Perspective','Callback', 'if strcmp(get(gca,''Projection''),''orthographic'')  set(gca,''Projection'',''perspective''); else set(gca,''Projection'',''orthographic''); end');
      uimenu(uicm, 'Separator','on','Label', 'About ResLibCal...', ...
        'Callback',[ 'msgbox(''' ResLibCal('version') ''',''About ResLibCal'',''help'')' ]);
      set(gca, 'UIContextMenu', uicm, 'Tag','ResLibCal_View3_Context');
    end
  end
