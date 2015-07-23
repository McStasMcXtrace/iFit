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

for index=1:numel(resolutions)
  resolution = resolutions{index};
  H=resolution.HKLE(1); K=resolution.HKLE(2); 
  L=resolution.HKLE(3); W=resolution.HKLE(4);

  if ~isempty(strfind(mode,'rlu'))
    NP       = resolution.abc.RM;
    FrameStr = resolution.abc.FrameStr;
    cloud    = resolution.abc.cloud;
    Labels   = {'Q_1','Q_2','Q_3','\omega'};
    Units    = 'rlu';
    centre   = resolution.abc.hkl2Frame*[ H K L ]'; centre(4) = W;
  else
    NP       = resolution.xyz.RM;
    FrameStr = resolution.xyz.FrameStr;
    cloud    = resolution.xyz.cloud;
    Labels   = {'Q_x','Q_y','Q_z','\omega'};
    Units    = 'Angs-1';
    centre   = resolution.xyz.hkl2Frame*[ H K L ]'; centre(4) = W;
  end
  FrameStr{end+1} = ''; % no specific axis for energy (4)

  if isempty(NP) || ~all(isreal(NP)), return; end
  
  if isempty(strfind(mode,'cloud')), cloud=[]; end
  
  % reduce dimensionality of the resolution mtrix
  if ~isempty(strfind(mode,'qz'))
    [dummy, NP] = rc_int(4,1, NP); % this function strips out the row=col=4, and corrects determinant
                                   % using the cofactor rule.
    ResLibCal_Proj_plot3D([1 2 3], NP, FrameStr, Labels, Units, cloud, centre);  % xyz
    title([ 'Resolution in ' Labels{1:3} ' - ' out.EXP.method ])
    if index < numel(resolutions), hold on; else hold off; end
  else
    [dummy, NP] = rc_int(3,1, NP);
    ResLibCal_Proj_plot3D([1 2 4], NP, FrameStr, Labels, Units, cloud, centre);  % xye
    title([ 'Resolution in ' Labels{[1 2 4]} ' - ' out.EXP.method ])
    if index < numel(resolutions), hold on; else hold off; end
  end
  
end % for

% ------------------------------------------------------------------------------
function h=ResLibCal_Proj_plot3D(index, NP, FrameStr, Labels, Units, cloud, centre)

  % plot the cloud projection if any
  if ~isempty(cloud)
    x=cloud{index(1)}; y=cloud{index(2)}; z=cloud{index(3)};
    if isempty(find(index == 4)), e = cloud{4};
    else e = cloud{3}; end
    if numel(x) > 200, x=x(1:200); end
    if numel(y) > 200, y=y(1:200); end
    if numel(z) > 200, z=z(1:200); end
    if numel(e) > 200, e=e(1:200); end
    % h=scatter3(x,y,z,3,e,'filled');
    h=plot3(x,y,z,'o');
    set(h,'DisplayName',[ Labels{index} ' (cloud)' ], 'Tag', 'ResLibCal_View3_Cloud');
    hold on
  end
  
  % plot the ellipse on top.
  [h, XX,YY,ZZ] = Ellipse_plot(NP, centre(index),10);
  set(h,'DisplayName',[ Labels{index} ],'Tag','ResLibCal_View3_Volume');
  
  sx=max(XX(:)) - min(XX(:)); ix=index(1);
  sy=max(YY(:)) - min(YY(:)); iy=index(2);
  sz=max(ZZ(:)) - min(ZZ(:)); iz=index(3);
  
  xlabel({[ '{\bf ' Labels{ix} '} ' FrameStr{ix} ' [' Units ']' ], ...
       [ '{\delta}' Labels{ix} '=' num2str(sx, 3) ]})
  ylabel({[ '{\bf ' Labels{iy} '} ' FrameStr{iy} ' [' Units ']' ], ...
       [ '{\delta}' Labels{iy} '=' num2str(sy, 3) ]})
  zlabel({[ '{\bf ' Labels{iz} '} ' FrameStr{iz} ' [' Units ']' ], ...
       [ '{\delta}' Labels{iz} '=' num2str(sy, 3) ]})
       
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
        'clear tmp_a tmp_e; lighting none;' ]);
    uimenu(uicm, 'Label','Add Light','Callback', 'light;lighting phong;');
    uimenu(uicm, 'Label','Transparency','Callback', 'alpha(0.5);');
    uimenu(uicm, 'Label','Toggle Perspective','Callback', 'if strcmp(get(gca,''Projection''),''orthographic'')  set(gca,''Projection'',''perspective''); else set(gca,''Projection'',''orthographic''); end');
    uimenu(uicm, 'Separator','on','Label', 'About ResLibCal...', ...
      'Callback',[ 'msgbox(''' ResLibCal('version') ''',''About ResLibCal'',''help'')' ]);
    set(gca, 'UIContextMenu', uicm, 'Tag','ResLibCal_View3_Context');
  end
