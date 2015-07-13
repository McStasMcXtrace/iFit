function out=rc_projs(out, mode)
%
% MATLAB function to plot the projections of the resolution ellipse
% of a triple axis
%
% Input:
%  out:  EXP ResLib structure 
%  mode: can be set to 'rlu' so that the plot is in lattice RLU frame
%
% DFM 10.11.95

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
  disp([ mfilename ': Input argument does not contain EXP structure. Skipping.' ])
  return
end

% handle resolution for scans
if  numel(out.resolution) > 1
  for index=1:numel(out.resolution)
    this = out; 
    if iscell(this.resolution), this.resolution=this.resolution{index};
    else                        this.resolution=this.resolution(index); 
    end
    if index > 1, 
      this_mode = [ mode ' scan' ]; hold on
    else 
      this_mode=mode; hold off; 
    end
    feval(mfilename, this, this_mode);
  end
  hold off
  return
end

if isfield(out.resolution,'HKLE')
  H=out.resolution.HKLE(1); K=out.resolution.HKLE(2); 
  L=out.resolution.HKLE(3); W=out.resolution.HKLE(4);
  EXP.QH=H; EXP.QK=K; EXP.QL=L; EXP.W=W; % update single scan step
else
  EXP.QH=EXP.QH(1); EXP.QK=EXP.QK(1); EXP.QL=EXP.QL(1); EXP.W=EXP.W(1); 
end

R0  = out.resolution.R0;

if ~isempty(strfind(mode,'rlu'))
  NP  = out.resolution.RMS;

  [xvec,yvec,zvec,sample,rsample]=StandardSystem(EXP);

  %find reciprocal-space directions of X and Y axes

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
  NP = out.resolution.RM;
  frame = '[Qx,Qy,';
end
if isempty(strfind(mode,'qz')), frame = [ frame 'E]' ];
else frame = [ frame 'Qz]' ]; end
if isempty(NP) || ~all(isreal(NP)), return; end
A=NP;
const=1.17741; % half width factor

%----- Remove the vertical component from the matrix.
B = A;
if isempty(strfind(mode,'qz'))
  B(3,:)=[]; B(:,3)=[]; % remove Qz component
else
  B(4,:)=[]; B(:,4)=[]; % remove E component
end
% B=[A(1,1:2),A(1,4);A(2,1:2),A(2,4);A(4,1:2),A(4,4)];

%----- Work out projections for different cuts through the ellipse

%----- S is the rotation matrix that diagonalises the projected ellipse

%----- 1. Qx, Qy plane

[R0P,MP]=rc_int(3,R0,B);

theta=0.5*atan2(2*MP(1,2),(MP(1,1)-MP(2,2)));
S=[cos(theta) sin(theta); -sin(theta) cos(theta)];

MP=S*MP*S';

hwhm_xp=const/sqrt(MP(1,1));
hwhm_yp=const/sqrt(MP(2,2));
[x,y]=rc_ellip(hwhm_xp,hwhm_yp,theta); pqx=max(x)-min(x); pqy=max(y)-min(y);
subplot(2,2,1);
set(line(x,y),'Color','k')
if isempty(strfind(mode, 'scan')), fill(x,y,'r'); end
if ~isempty(strfind(mode,'rlu'))
  xlabel({[ '{\bf Q_1} (' o1 ') [rlu]' ], ...
    [ '{\delta}Q_1=' num2str(max(x)-min(x), 3) ]})
  ylabel({[ '{\bf Q_2} (' o2 ') [rlu]' ], ...
    [ '{\delta}Q_2=' num2str(max(y)-min(y), 3) ]})
else
  xlabel([ '{\bf Q_x} (long.) [A^{-1}] {\delta}Q_x=' num2str(max(x)-min(x),3) ] )
  ylabel([ '{\bf Q_y} (trans.) [A^{-1}] {\delta}Q_y=' num2str(max(y)-min(y),3) ] )
end
title(EXP.method); 
da=daspect; da(1:2) = max(da(1:2)); daspect(da);
pb=pbaspect; pb(1:2)=da(1); pbaspect(pb);
x1=xlim; x2=ylim;
grid on

%---------------- Add slice through Qx,Qy plane ----------------------

MP=A(1:2,1:2);

theta=0.5*atan2(2*MP(1,2),(MP(1,1)-MP(2,2)));
S=[cos(theta) sin(theta); -sin(theta) cos(theta)];

MP=S*MP*S';

hwhm_xp=const/sqrt(MP(1,1));
hwhm_yp=const/sqrt(MP(2,2));
[x,y]=rc_ellip(hwhm_xp,hwhm_yp,theta); dqx = max(x)-min(x); dqy=max(y)-min(y);
set(line(x,y),'LineStyle','--')

%----- 2. Qx, W plane

[R0P,MP]=rc_int(2,R0,B);

theta=0.5*atan2(2*MP(1,2),(MP(1,1)-MP(2,2)));
S=[cos(theta) sin(theta); -sin(theta) cos(theta)];

MP=S*MP*S';

hwhm_xp=const/sqrt(MP(1,1));
hwhm_yp=const/sqrt(MP(2,2)); 
[x,y]=rc_ellip(hwhm_xp,hwhm_yp,theta); pe=max(y)-min(y);
subplot(2,2,2);
set(line(x,y),'Color','k')
if isempty(strfind(mode, 'scan')), fill(x,y,'r'); end
if ~isempty(strfind(mode,'rlu'))
  xlabel({[ '{\bf Q_1} (' o1 ') [rlu]' ], ...
    [ '{\delta}Q_1=' num2str(max(x)-min(x),3) ]})
else
  xlabel([ '{\bf Q_x} (long.) [A^{-1}] {\delta}Q_x=' num2str(max(x)-min(x),3) ])
end

if isempty(strfind(mode,'qz'))
  ylabel([ '{\bf \omega} [meV]  {\delta}\omega=' num2str(max(y)-min(y),2) ])
else
  if ~isempty(strfind(mode,'rlu'))
    ylabel([ '{\bf Qz} (' o3 ') [rlu] - {\delta}Qz=' num2str(max(y)-min(y),3) ]);
  else
    ylabel([ '{\bf Qz} (vert.) [A^{-1}] - {\delta}Qz=' num2str(max(y)-min(y),3) ]);
  end
end
xlim(x1); xe=ylim;
grid on

%---------------- Add slice through Qx,W plane ----------------------

if isempty(strfind(mode,'qz'))
  MP=[A(1,1) A(1,4);A(4,1) A(4,4)];
else
  MP=[A(1,1) A(1,3);A(3,1) A(3,3)];
end

theta=0.5*atan2(2*MP(1,2),(MP(1,1)-MP(2,2)));
S=[cos(theta) sin(theta); -sin(theta) cos(theta)];

MP=S*MP*S';

hwhm_xp=const/sqrt(MP(1,1));
hwhm_yp=const/sqrt(MP(2,2));
[x,y]=rc_ellip(hwhm_xp,hwhm_yp,theta); de=max(y)-min(y);
set(line(x,y),'LineStyle','--')

if isempty(strfind(mode,'qz'))
  title('Energy resolution')
else
  title('Vertical Qz resolution')
end

%----- 3. Qy, W plane

[R0P,MP]=rc_int(1,R0,B);

theta=0.5*atan2(2*MP(1,2),(MP(1,1)-MP(2,2)));
S=[cos(theta) sin(theta); -sin(theta) cos(theta)];

MP=S*MP*S';

hwhm_xp=const/sqrt(MP(1,1));
hwhm_yp=const/sqrt(MP(2,2));
[x,y]=rc_ellip(hwhm_xp,hwhm_yp,theta);
subplot(2,2,3);
set(line(x,y),'Color','k')
if isempty(strfind(mode, 'scan')), fill(x,y,'r'); end
if ~isempty(strfind(mode,'rlu'))
  xlabel({[ '{\bf Q_2} (' o2 ') [rlu]' ], ...
    [ '{\delta}Q_2=' num2str(max(x)-min(x),3) ]})
else
  xlabel([ '{\bf Q_y} (trans.) [A^{-1}] {\delta}Q_y=' num2str(max(x)-min(x)) ])
end

if isempty(strfind(mode,'qz'))
  ylabel([ '{\bf \omega} [meV]  {\delta}\omega=' num2str(max(y)-min(y),2) ])
else
  if ~isempty(strfind(mode,'rlu'))
    ylabel([ '{\bf Qz} (' o3 ') [rlu] - {\delta}Qz=' num2str(max(y)-min(y),3) ]);
  else
    ylabel([ '{\bf Qz} (vert.) [A^{-1}] - {\delta}Qz=' num2str(max(y)-min(y),3) ]);
  end
end

xlim(x2); ylim(xe);
title([ 'ResLibCal ' datestr(now) ])
grid on

%---------------- Add slice through Qy,W plane ----------------------
if isempty(strfind(mode,'qz'))
  MP=[A(2,2) A(2,4);A(4,2) A(4,4)];
else
  MP=[A(2,2) A(2,3);A(3,2) A(3,3)];
end

theta=0.5*atan2(2*MP(1,2),(MP(1,1)-MP(2,2)));
S=[cos(theta) sin(theta); -sin(theta) cos(theta)];

MP=S*MP*S';

hwhm_xp=const/sqrt(MP(1,1));
hwhm_yp=const/sqrt(MP(2,2));
[x,y]=rc_ellip(hwhm_xp,hwhm_yp,theta);
set(line(x,y),'LineStyle','--')

% display a text edit uicontrol so that users can select/copy/paste

[res, inst] = ResLibCal_FormatString(out, mode);

message = [ res, inst ];

% fill 4th sub-panel with uicontrol
p(1) = 0.5; p(2) = 0.01; p(3) = 0.49; p(4) = 0.49;
h = findobj(gcf, 'Tag','ResLibCal_View2_Edit');
if isempty(h)
  h = uicontrol('Tag','ResLibCal_View2_Edit', ...
    'Style','edit','Units','normalized',...
    'Position',p, 'Max',2,'Min',0, ...
    'String', message, 'FontName', 'FixedWidth', ...
    'FontSize',8,'BackgroundColor','white','HorizontalAlignment','left');
else
  set(h, 'String',message);
end

hold off

% ==============================================================================

function  [x,y] = rc_ellip(a,b,phi,x0,y0,n)
% ELLIPSE  Plotting ellipse.
%       ELLIPSE(A,B,PHI,X0,Y0,N)  Plots ellipse with
%	semiaxes A, B, rotated by the angle PHI,
%	with origin at X0, Y0 and consisting of N points
%	(default 100).
%	[X,Y] = ELLIPSE(...) Instead of plotting returns
%	coordinates of the ellipse.

%  Kirill K. Pankratov, kirill@plume.mit.edu
%  03/21/95

n_dflt = 36;

if nargin < 6, n = n_dflt; end
if nargin < 5, y0 = 0; end
if nargin < 4, x0 = 0; end
if nargin < 3, phi = 0; end
if nargin < 2, b = 1; end
if nargin < 1, a = 1; end


th = linspace(0,2*pi,n+1);
x = a*cos(th);
y = b*sin(th);

c = cos(phi);
s = sin(phi);

th = x*c-y*s+x0;
y = x*s+y*c+y0;
x = th;

if nargout==0, plot(x,y); end

