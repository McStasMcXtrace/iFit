function filename=iData_private_saveas_vtk(a, filename)
% private function to write VTK files

% inline: WriteToVTK, export3Dline2VTK, writeVTK

if ndims(a) ~= 2 & ndims(a) ~= 3
  filename=[];
  iData_private_warning(mfilename,[ 'Can only export 2D and 3D objects to VTK format. Object ' a.Tag ' has ndims=' num2str(ndims(a)) ]);
  return
end

NL = sprintf('\n');
str = [ '# URL: ifit.mccode.org' NL ...
        '# Creator: iFit/@iData/saveas - ' version(a) ];

m = get(a,'Monitor'); 
if not(all(m==0| m==1))
  title(a, [ title(a) ' per monitor' ]);
  m=mean(m(:));
  str = [ str NL '# SignalNormalizedToMonitor: Yes' NL ...
                 '# MeanMonitor: ' num2str(m) NL ];
end

x = getaxis(a, 2);
y = getaxis(a, 1);

str = [ str NL '# Title: ' a.Title NL ...
          '# Label: ' a.Label NL ...
          '# DisplayName: ' a.DisplayName NL ...
          '# User: ' a.User NL ...
          '# CreationDate: ' a.Date NL ...
          '# ModificationDate: ' a.ModificationDate NL ...
          '# Tag: ' a.Tag NL ...
          '# Axis2_Label' xlabel(a) NL ...
          '# Axis2_Min' num2str(min(x(:))) NL ...
          '# Axis2_Max', num2str(max(x(:))) NL ...
          '# Axis2_Length', num2str(length(unique(x(:)))) NL ...
          '# Axis1_Label' ylabel(a) NL ...
          '# Axis1_Min', num2str(min(y(:))) NL ...
          '# Axis1_Max', num2str(max(y(:))) NL ...
          '# Axis1_Length', num2str(length(unique(y(:)))) NL ];
          
str=[ char(a) NL ];

if ndims(a) == 3 && ~isvector(a)
  Signal = getaxis(a, 0); % Signal/Monitor
  WriteToVTK(Signal, filename, str); % dump the matrix directly to the file as POINTS
elseif ndims(a) == 2 && isvector(a)
  export3Dline2VTK(filename, [ x(:) y(:) Signal(:) ], str); % export as VTK3 LINES, can not re-import
elseif ndims(a) == 2
  a = interp(a,'grid');
  x = getaxis(a, 2);
  y = getaxis(a, 1);
  x =[reshape(x,numel(x),1),reshape(y,numel(y),1)];
  writeVTK(filename, x,delaunayn(x),getaxis(a,0), str);
elseif ndims(a) == 3
  a = interp(a,'grid');
  x = getaxis(a, 2);
  y = getaxis(a, 1);
  z = getaxis(a, 3);
  p=[ x(:) y(:) z(:) ]; 
  u = getaxis(a,0); u=u(:);
  writeVTK(filename, p,delaunayn(p),u, str);
end

% ------------------------------------------------------------------------------

function WriteToVTK(matrix, filename, str)
% WriteToVTK(matrix, filename)
%
% Writes a 3D matrix as a VTK file. View with paraview.
% The matrix must be 3D and is normalised (kind of).
%
% A sligthly old format is used because it is simpler.

% Example VTK file:
% # vtk DataFile Version 2.0
% Volume example
% ASCII
% DATASET STRUCTURED_POINTS
% DIMENSIONS 3 4 6
% ASPECT_RATIO 1 1 1
% ORIGIN 0 0 0
% POINT_DATA 72
% SCALARS volume_scalars char 1
% LOOKUP_TABLE default
% 0 0 0 0 0 0 0 0 0 0 0 0
% 0 5 10 15 20 25 25 20 15 10 5 0
% 0 10 20 30 40 50 50 40 30 20 10 0
% 0 10 20 30 40 50 50 40 30 20 10 0
% 0 5 10 15 20 25 25 20 15 10 5 0
% 0 0 0 0 0 0 0 0 0 0 0 0

% Get the matrix dimensions.
[N M O] = size(matrix);

% Get the maximum value for the normalisation.
mx = max(matrix(:));

% Open the file.
fid = fopen(filename, 'w');
if fid == -1
    error('Cannot open file for writing.');
end

% New line.
nl = sprintf('\n'); % Stupid matlab doesn't interpret \n normally.

% Write the file header.
fwrite(fid, ['# vtk DataFile Version 2.0' nl str 'ASCII' nl ...
    'DATASET STRUCTURED_POINTS' nl 'DIMENSIONS ' ...
    num2str(N) ' ' num2str(M) ' ' num2str(O) nl 'ASPECT_RATIO 1 1 1' nl ...
    'ORIGIN 0 0 0' nl 'POINT_DATA ' ...
    num2str(N*M*O) nl 'SCALARS volume_scalars char 1' nl 'LOOKUP_TABLE default' nl]);

for z = 1:O
    % Get this layer.
    v = matrix(:, :, z);
    % Scale it. This assumes there are no negative numbers. I'm not sure
    % this is actually necessary.
    v = round(100 .* v(:) ./ mx);
    % Write the values as text numbers.
    fwrite(fid, num2str(v'));
    % Newline.
    fwrite(fid, nl);
    
    % Display progress.
    disp([num2str(round(100*z/O)) '%']);
end

% Close the file.
fclose(fid);

% end WriteToVTK

% ------------------------------------------------------------------------------

function export3Dline2VTK(file,path3D,str)
%
% This function takes as input a 3D line and export
% it to an ASCII VTK file which can be oppened with the viewer Paraview.
%
% Input :
%           "file" is the name without extension of the file (string).
%           "path3D" is a list of 3D points (nx3 matrix)
%           "dir" is the path of the directory where the file is saved (string). (Optional)
%
% Simple example :
%
%   Xi=linspace(0,10,200)';
%   Yi=[zeros(100,1);linspace(0,20,100)'];
%   Zi=5*sin(5*Xi);
%   export3Dline2VTK('ExampleLine',[Xi Yi Zi])
%
% David Gingras, January 2009

X=path3D(:,1);
Y=path3D(:,2);
Z=path3D(:,3);
nbWaypoint=length(X);
if mod(nbWaypoint,3)==1
    X(nbWaypoint+1:nbWaypoint+2,1)=[0;0];
    Y(nbWaypoint+1:nbWaypoint+2,1)=[0;0];
    Z(nbWaypoint+1:nbWaypoint+2,1)=[0;0];
elseif mod(nbWaypoint,3)==2
    X(nbWaypoint+1,1)=0;
    Y(nbWaypoint+1,1)=0;
    Z(nbWaypoint+1,1)=0;
end
nbpoint=length(X);

fid=fopen(file,'wt');
fprintf(fid,'# vtk DataFile Version 3.0\n%sASCII\nDATASET POLYDATA\n', str);
fprintf(fid,'POINTS %d float\n',nbpoint);
fprintf(fid,'%3.7f %3.7f %3.7f %3.7f %3.7f %3.7f %3.7f %3.7f %3.7f\n',[X(1:3:end-2) Y(1:3:end-2) Z(1:3:end-2) X(2:3:end-1) Y(2:3:end-1) Z(2:3:end-1) X(3:3:end) Y(3:3:end) Z(3:3:end)]');

if mod(nbWaypoint,2)==0
    nbLine=2*nbWaypoint-2;
else
    nbLine=2*(nbWaypoint-1);
end

ass1=zeros(nbLine,1);
ass2=zeros(nbLine,1);

ass1(1)=0;
ass2(1)=1;
ass1(end)=nbLine/2;
ass2(end)=nbLine/2-1;

ass1(2:2:nbLine-1)=1:nbLine/2-1;
ass1(3:2:nbLine)=1:nbLine/2-1;
ass2(2:2:nbLine-1)=ass1(2:2:nbLine-1)-1;
ass2(3:2:nbLine)=ass1(3:2:nbLine)+1;

fprintf(fid,'\nLINES %d %d\n',nbLine,3*nbLine);
fprintf(fid,'2 %d %d\n',ass1',ass2');

fclose(fid);

% end export3Dline2VTK

% ------------------------------------------------------------------------------

function writeVTK(filename,t,p,u, str)
% vtk export
% creates a vtk-file filename.vtk containing a simplicial mesh (2- or 3d)
% and additional vertex data
%
% input: filename   destination 
%                   (string)
%        p          array of N points 
%                   (Nxd matrix where d denotes the dimension)
%        t          triangulation/tetrahedralization of the points in p
%                   (Mxd+1 array, where M denotes the number of simplices)
%        u          mesh function assigning a real number to every point in
%                   p (the vertices of the
%                   triangulation/tetrahedralization)
%                   (Nx1 array)
%
% example usage:
%        2d: [X,Y]=meshgrid(0:0.01:1,0:0.01:1);
%            p=[reshape(X,prod(size(X)),1),reshape(Y,prod(size(Y)),1)];
%            t=delaunayn(p);
%            u=peaks(6*p(:,1)-3,6*p(:,2)-3);
%            writeVTK('test2d',t,p,u);
%        3d: p=rand(10,3); 
%            t=delaunayn(p); 
%            u=sum(p.^2,2);
%            writeVTKcell('test3d',t,p,u);
% (the result is accessible with paraview!)
%
% (c) Daniel Peterseim, 2009-11-07

[np,dim]=size(p);
[nt]=size(t,1);
celltype=[3,5,10];

FID = fopen(filename,'w+');
fprintf(FID,'# vtk DataFile Version 2.0\n%sASCII\n', str);
fprintf(FID,'DATASET UNSTRUCTURED_GRID\n');

fprintf(FID,'POINTS %d float\n',np);
s='%f %f %f \n';
P=[p zeros(np,3-dim)];
fprintf(FID,s,P');

fprintf(FID,'CELLS %d %d\n',nt,nt*(dim+2));
s='%d ';
for k=1:dim+1
    s=horzcat(s,{' %d'});
end
s=cell2mat(horzcat(s,{' \n'}));
fprintf(FID,s,[(dim+1)*ones(nt,1) t-1]');

fprintf(FID,'CELL_TYPES %d\n',nt);
s='%d\n';
fprintf(FID,s,celltype(dim)*ones(nt,1));

fprintf(FID,'POINT_DATA %s\nSCALARS u float 1\nLOOKUP_TABLE default\n',num2str(np));
s='%f\n';
fprintf(FID,s,u);

fclose(FID);

% end writeVTK
