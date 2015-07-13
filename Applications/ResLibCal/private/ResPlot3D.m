function out = ResPlot3D(out, mode)
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
  delete(findobj(gcf,'Tag','ResLibCal_View3_Proj1'));
  delete(findobj(gcf,'Tag','ResLibCal_View3_Proj2'));
  delete(findobj(gcf,'Tag','ResLibCal_View3_Proj3'));
  delete(findobj(gcf,'Tag','ResLibCal_View3_Volume'));
  delete(findobj(gcf,'Tag','ResLibCal_View3_Cloud'));
end

% handle resolution for scans
if  numel(out.resolution) > 1
  for index=1:numel(out.resolution)
    this = out; 
    if iscell(this.resolution), this.resolution=this.resolution{index};
    else                        this.resolution=this.resolution(index);
    end
    if index > 1, this_mode = [ mode ' scan' ]; else this_mode=mode; end
    feval(mfilename, this, this_mode);
    hold on
  end
  hold off
  return
end

% handle single scan step
if isfield(out.resolution,'HKLE')
  H=out.resolution.HKLE(1); K=out.resolution.HKLE(2); 
  L=out.resolution.HKLE(3); W=out.resolution.HKLE(4);
else
  H=EXP.QH(1); K=EXP.QK(1); L=EXP.QL(1); W=EXP.W(1);
end
EXP.QH=H; EXP.QK=K; EXP.QL=L; EXP.W=W; % update single scan step

if ~isempty(strfind(mode,'rlu'))
  RMS  = out.resolution.RMS;
  [xvec,yvec,zvec,sample,rsample]=StandardSystem(EXP);
  qx=scalar(xvec(1,:),xvec(2,:),xvec(3,:),H,K,L,rsample);
  qy=scalar(yvec(1,:),yvec(2,:),yvec(3,:),H,K,L,rsample);

  % Q vectors on figure axes
  o1=EXP.orient1;
  o2=EXP.orient2;
  pr=scalar(o2(1),o2(2),o2(3),yvec(1),yvec(2),yvec(3),rsample);
  o2 = yvec*pr;

  o1 = o1(:)';
  o2 = o2(:)';
  o3 = cross(o1,o2);
  
  % convert o1 and o2 to normalised strings
  if all(abs(o1 - round(o1)) < 1e-5)
    o1 = round(o1);
    if sum(o1.*o1) > 1,  o1=[ '[' sprintf('%i ', o1) ']/\surd' num2str(sum(o1.*o1)) ];
    else                 o1=[ '[' sprintf('%i ', o1) ']' ];
    end
  else o1 = [ '[' sprintf('%.3f ', o1) ']' ]; end

  if all(abs(o2 - round(o2)) < 1e-5)
    o2 = round(o2);
    if sum(o2.*o2) > 1,  o2=[ '[' sprintf('%i ', o2) ']/\surd' num2str(sum(o2.*o2)) ];
    else                 o2=[ '[' sprintf('%i ', o2) ']' ];
    end
  else o2 = [ '[' sprintf('%.3f ', o2) ']' ]; end
  
  if all(abs(o3 - round(o3)) < 1e-5)
    o3 = round(o3);
    if sum(o3.*o3) > 1,  o3=[ '[' sprintf('%i ', o3) ']/\surd' num2str(sum(o3.*o3)) ];
    else                 o3=[ '[' sprintf('%i ', o3) ']' ];
    end
  else o3 = [ '[' sprintf('%.3f ', o3) ']' ]; end
  
  frame = '[Q1,Q2,';
else
  qx = 0; qy=0;
  RMS  = out.resolution.RM;
  
  frame = '[Qx,Qy,';
end
if isempty(strfind(mode,'qz')),
  frame = [ frame 'E]' ];
  qw=W;
else
  frame = [ frame 'Qz]' ];
  qw=L;
end
if isempty(RMS) || ~all(isreal(RMS)), return; end
SMAGridPoints      =40;
EllipsoidGridPoints=40;

len=1;
center=1;

% for this plot to work, we need to remove row-column 3 of RMS
if isempty(strfind(mode,'qz')) % plot(Q1,Q2,W), remove Qz=3
  RMS(:,3) = [];
  RMS(3,:) = [];
else                           % plot(Q1,Q2,Qz), remove E=4
  RMS(:,4) = [];
  RMS(4,:) = [];
end

hold on

%plot ellipsoids
wx=fproject(RMS,1); wy=fproject(RMS,2); ww=fproject(RMS,3);
for point=1:len
   x=linspace(-wx(point)*1.5,wx(point)*1.5,EllipsoidGridPoints)+qx(point);
   y=linspace(-wy(point)*1.5,wy(point)*1.5,EllipsoidGridPoints)+qy(point);
   z=linspace(-ww(point)*1.5,ww(point)*1.5,EllipsoidGridPoints)+qw(point);
   [xg,yg,zg]=meshgrid(x,y,z);
   ee=RMS(1,1,point)*(xg-qx(point)).^2+RMS(2,2,point)*(yg-qy(point)).^2+RMS(3,3,point)*(zg-qw(point)).^2+...
      2*RMS(1,2,point)*(xg-qx(point)).*(yg-qy(point))+...
      2*RMS(1,3,point)*(xg-qx(point)).*(zg-qw(point))+...
      2*RMS(3,2,point)*(zg-qw(point)).*(yg-qy(point));
   p = patch(isosurface(xg, yg, zg, ee, 2*log(2)));
   isonormals(xg,yg,zg,ee,p)
   set(p, 'FaceColor', 'red', ...
          'EdgeColor', 'none','BackFaceLighting','reverselit', ...
          'Tag','ResLibCal_View3_Volume','DisplayName','ellipsoid');
end;

%da=daspect;
%daspect([da(2) da(2) da(3)]);
RANGE = [ xlim ylim zlim ];

%plot projections (not for scans as it makes plot crowded)
if isempty(strfind(mode,'scan'))
  [proj3,sec]=project(RMS,3);
  [proj2,sec]=project(RMS,2);
  [proj1,sec]=project(RMS,1);
  phi=0.1:2*pi/3000:2*pi+0.1;
  for i=1:len
     r3=sqrt(2*log(2)./(proj3(1,1,i)*cos(phi).^2+proj3(2,2,i)*sin(phi).^2+2*proj3(1,2,i)*cos(phi).*sin(phi)));
     r2=sqrt(2*log(2)./(proj2(1,1,i)*cos(phi).^2+proj2(2,2,i)*sin(phi).^2+2*proj2(1,2,i)*cos(phi).*sin(phi)));
     r1=sqrt(2*log(2)./(proj1(1,1,i)*cos(phi).^2+proj1(2,2,i)*sin(phi).^2+2*proj1(1,2,i)*cos(phi).*sin(phi)));
     xproj3=r3.*cos(phi)+qx(i);   yproj3=r3.*sin(phi)+qy(i);   zproj3=ones(size(xproj3))*RANGE(5);
     xproj2=r2.*cos(phi)+qx(i);   zproj2=r2.*sin(phi)+qw(i);   yproj2=ones(size(xproj2))*RANGE(4);
     yproj1=r1.*cos(phi)+qy(i);   zproj1=r1.*sin(phi)+qw(i);   xproj1=ones(size(yproj1))*RANGE(2);
     h=plot3(xproj1,yproj1,zproj1); 
     set(h,'Tag','ResLibCal_View3_Proj1','DisplayName','proj_x');
     h=plot3(xproj2,yproj2,zproj2); 
     set(h,'Tag','ResLibCal_View3_Proj2','DisplayName','proj_y');
     h=plot3(xproj3,yproj3,zproj3); 
     set(h,'Tag','ResLibCal_View3_Proj3','DisplayName','proj_z');
  end;
end

da=daspect; da(1:2) = max(da(1:2)); daspect(da);

% add contextual menu
if isempty(findobj(gcf,'Tag','ResLibCal_View3_Context'))
  %finalize 3D plot
  box on; grid on;
  view(3);

  uicm = uicontextmenu;
  uimenu(uicm, 'Label','Toggle grid', 'Callback','grid');
  uimenu(uicm, 'Label','Reset Flat/3D View', 'Callback', [ ...
      '[tmp_a,tmp_e]=view; if (tmp_a==0 & tmp_e==90) view(3); else view(2); end;' ...
      'clear tmp_a tmp_e; lighting none;' ]);
    uimenu(uicm, 'Label','Add Light','Callback', 'light;lighting phong;');
    uimenu(uicm, 'Label','Transparency','Callback', 'alpha(0.7);');
  set(gca, 'UIContextMenu', uicm, 'Tag','ResLibCal_View3_Context');
    uimenu(uicm, 'Label','Toggle Perspective','Callback', 'if strcmp(get(gca,''Projection''),''orthographic'')  set(gca,''Projection'',''perspective''); else set(gca,''Projection'',''orthographic''); end');
end
% add labels
if ~isempty(strfind(mode,'rlu'))
  xlabel({[ '{\bf Q_1} (' o1 ') [rlu]'],[ '{\delta}Q_1=' num2str(max(x)-min(x)) ]})
  ylabel({[ '{\bf Q_2} (' o2 ') [rlu]'],[ '{\delta}Q_2=' num2str(max(y)-min(y)) ]})
else
  xlabel([ '{\bf Q_x} [A^{-1}] {\delta}Q_x=' num2str(max(x)-min(x)) ])
  ylabel([ '{\bf Q_y} [A^{-1}] {\delta}Q_y=' num2str(max(y)-min(y)) ])
end
if isempty(strfind(mode,'qz'))
  zlabel([ '{\bf \omega} (meV) - {\delta}\omega=' num2str(max(z)-min(z)) ]);
else
  if ~isempty(strfind(mode,'rlu'))
    zlabel([ '{\bf Qz} (' o3 ') [rlu] - {\delta}Qz=' num2str(max(z)-min(z)) ]);
  else
    zlabel([ '{\bf Qz} [A^{-1}] - {\delta}Qz=' num2str(max(z)-min(z)) ]);
  end
end
title([ 'Resolution in ' frame ' - ' out.EXP.method ])

% plot cloud of points if available
if ~isempty(strfind(mode,'rlu')) && isfield(out.resolution, 'cloud')
  if isempty(strfind(mode,'qz')), ax=[1 2 4]; else ax=[1 2 3]; end
  h=plot3(out.resolution.cloud{ax},'.');
  set(h,'Tag','ResLibCal_View3_Cloud','DisplayName','cloud','MarkerSize',0.5);
end


hold off

%========================================================================================================
%========================================================================================================
%========================================================================================================

function hwhm=fproject (mat,i)



if (i==1) v=3;j=2;end;
if (i==2) v=1;j=3;end;
if (i==3) v=1;j=2;end;
[a,b,c]=size(mat);
proj=zeros(2,2,c);
proj(1,1,:)=mat(i,i,:)-mat(i,v,:).^2./mat(v,v,:);
proj(1,2,:)=mat(i,j,:)-mat(i,v,:).*mat(j,v,:)./mat(v,v,:);
proj(2,1,:)=mat(j,i,:)-mat(j,v,:).*mat(i,v,:)./mat(v,v,:);
proj(2,2,:)=mat(j,j,:)-mat(j,v,:).^2./mat(v,v,:);

hwhm=proj(1,1,:)-proj(1,2,:).^2./proj(2,2,:);
hwhm=sqrt(2*log(2))./sqrt(hwhm);


function [proj, sec]=project (mat,v)
if v == 3, i=1;j=2; end;
if v == 1, i=2;j=3; end;
if v == 2, i=1;j=3; end;
[a,b,c]=size(mat);
proj=zeros(2,2,c);
sec=zeros(2,2,c);
proj(1,1,:)=mat(i,i,:)-mat(i,v,:).^2./mat(v,v,:);
proj(1,2,:)=mat(i,j,:)-mat(i,v,:).*mat(j,v,:)./mat(v,v,:);
proj(2,1,:)=mat(j,i,:)-mat(j,v,:).*mat(i,v,:)./mat(v,v,:);
proj(2,2,:)=mat(j,j,:)-mat(j,v,:).^2./mat(v,v,:);
sec(1,1,:)=mat(i,i,:);
sec(1,2,:)=mat(i,j,:);
sec(2,1,:)=mat(j,i,:);
sec(2,2,:)=mat(j,j,:);
