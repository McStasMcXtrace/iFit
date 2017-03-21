function filename=iData_private_saveas_vtk(a, filename)
% private function to write VTK files

% inline: WriteToVTK (3D matrix, no axes, not used but keep it here), 
% export3Dline2VTK (xyz line), writeVTK (3D object, binary), exportTriangulation2VTK (2D surface), WriteToVTK(3D with axes, ascii)

if ndims(a) ~= 2 & ndims(a) ~= 3
  filename=[];
  warning([ mfilename, ': Can only export 2D and 3D objects to VTK format.\n\tObject ' a.Tag ' has ndims=' num2str(ndims(a)) ]);
  return
end

a = iData_private_cleannaninf(a);
str=char(a);

if ndims(a) <= 2 && isvector(a) % line, 1D or 2D
  if ndims(a) == 2
    x = getaxis(a, 2);
    y = getaxis(a, 1);
  else
    x = getaxis(a, 1);
    y = zeros(size(x));
  end
  Signal = getaxis(a, 0);
  export3Dline2VTK(filename, [ x(:) y(:) Signal(:) ], str); % export as VTK3 LINES
elseif ndims(a) == 2 % surface with axes
  a = interp(a,'grid');
  x = getaxis(a, 2); x=x(:);
  y = getaxis(a, 1); y=y(:);
  z = getaxis(a, 0); z=z(:);
  exportTriangulation2VTK(filename, [x y z], delaunay(x,y), str);
elseif ndims(a) == 3 % volume (no axes)
  a = interp(a,'grid');
  try
    vtkwrite(filename, 'structured_grid', ...
     getaxis(a, 1), getaxis(a, 2), getaxis(a, 3), ...
     'scalars','iData', getaxis(a, 0));
  catch
    WriteToVTK(getaxis(a,0), filename, str);
  end
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

if M*N*O > 1e5 % if big data set, we save in binary
  fmt = 'BINARY';
else
  fmt = 'ASCII';
end

% Get the maximum value for the normalisation.
mx = max(matrix(:));

% Open the file.
fid = fopen(filename, 'w');
if fid == -1
    error([ 'WriteToVTK: Cannot open file ' filename ' for writing.' ]);
end

% New line.
nl = sprintf('\n'); % Stupid matlab doesn't interpret \n normally.

% Write the file header.
fwrite(fid, ['# vtk DataFile Version 2.0' nl ...
    str nl ...
    fmt nl ...
    'DATASET STRUCTURED_POINTS' nl ...
    'DIMENSIONS ' num2str(N) ' ' num2str(M) ' ' num2str(O) nl ...
    'ASPECT_RATIO 1 1 1' nl ...
    'ORIGIN 0 0 0' nl ....
    'POINT_DATA '  num2str(N*M*O) nl ...
    'SCALARS volume_scalars char 1' nl ] );

if strcmp(fmt, 'ASCII')
  fwrite(fid, [ 'LOOKUP_TABLE default' nl ]);
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
else % BINARY
  tp = class(matrix);
  if( strcmp(tp, 'uint8') > 0 )      typ = 'unsigned_char';
  elseif( strcmp(tp, 'int8') > 0 )   typ = 'char';
  elseif( strcmp(tp, 'uint16') > 0 ) typ = 'unsigned_short';
  elseif( strcmp(tp, 'int16') > 0 )  typ = 'short';
  elseif( strcmp(tp, 'uint32') > 0 ) typ = 'unsigned_int';
  elseif( strcmp(tp, 'int32') > 0 )  typ = 'int';
  elseif( strcmp(tp, 'single') > 0 ) typ = 'float';
  elseif( strcmp(tp, 'double') > 0 ) typ = 'double';
  end
  fprintf(fid, 'SCALARS image_data %s\n', typ);
  fprintf(fid, 'LOOKUP_TABLE default\n');
  fwrite(fid,matrix,class(matrix));      % write data as binary
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

%  Copyright (c) 2009, David Gingras
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are
%  met:
%
%      * Redistributions of source code must retain the above copyright
%        notice, this list of conditions and the following disclaimer.
%      * Redistributions in binary form must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in
%        the documentation and/or other materials provided with the distribution
%
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%  POSSIBILITY OF SUCH DAMAGE.

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
fprintf(fid,'# vtk DataFile Version 3.0\n%s\nASCII\nDATASET POLYDATA\n', str);
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

function exportTriangulation2VTK(file,XYZ,tri,str)
%
% This function takes as input a 2D unrestricted triangulation and export
% it to an ASCII VTK file which can be oppened with the viewer Paraview.
%
% Input :
%           "dir" is the path of the directory where the file is saved (string). (Optional)
%           "file" is the name without extension of the file (string).
%           "XYZ" is the coordinate of the vertex of the triangulation (nx3 matrix).
%           "tri" is the list of triangles which contain indexes of XYZ (mx3 matrix).
%
% Sample example :
%
%   [X,Y,Z]=peaks(25);
%   X=reshape(X,[],1);
%   Y=reshape(Y,[],1);
%   Z=0.4*reshape(Z,[],1);
%   tri = delaunay(X,Y);
%   exportTriangulation2VTK('sampleExampleTri',[X Y Z],tri)
%
% Note : If the triangulation doesn't have Z component (a plane), put the
% third column of XYZ with all zeros. Paraview only deals with 3D object.
%
% David Gingras, January 2009

%  Copyright (c) 2009, David Gingras
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are
%  met:
%
%      * Redistributions of source code must retain the above copyright
%        notice, this list of conditions and the following disclaimer.
%      * Redistributions in binary form must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in
%        the documentation and/or other materials provided with the distribution
%
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%  POSSIBILITY OF SUCH DAMAGE.

X=XYZ(:,1);
Y=XYZ(:,2);
Z=XYZ(:,3);
nbpoint=length(X);
if mod(nbpoint,3)==1
    X(nbpoint+1:nbpoint+2,1)=[0;0];
    Y(nbpoint+1:nbpoint+2,1)=[0;0];
    Z(nbpoint+1:nbpoint+2,1)=[0;0];
elseif mod(nbpoint,3)==2
    X(nbpoint+1,1)=0;
    Y(nbpoint+1,1)=0;
    Z(nbpoint+1,1)=0;
end
nbpoint=length(X);
fid=fopen(file,'wt');
fprintf(fid,'# vtk DataFile Version 3.0\n%s\nASCII\nDATASET POLYDATA\n',str);
fprintf(fid,'POINTS %d float\n',nbpoint);
fprintf(fid,'%3.7f %3.7f %3.7f %3.7f %3.7f %3.7f %3.7f %3.7f %3.7f\n',[X(1:3:end-2) Y(1:3:end-2) Z(1:3:end-2) X(2:3:end-1) Y(2:3:end-1) Z(2:3:end-1) X(3:3:end) Y(3:3:end) Z(3:3:end)]');
ntri=length(tri);
fprintf(fid,'POLYGONS %d %d\n',ntri,4*ntri);

fprintf(fid,'3 %d %d %d\n',(tri-ones(ntri,3))');
fclose(fid);

% ------------------------------------------------------------------------------
%Copyright (c) 2015, Joe Yeh
%Copyright (c) 2009, David Gingras
%Copyright (c) 2011, W Thielicke
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are
%met:
%
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in
%      the documentation and/or other materials provided with the distribution
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%POSSIBILITY OF SUCH DAMAGE.
%
% <http://www.mathworks.com/matlabcentral/fileexchange/47814-export-3d-data-to-paraview-in-vtk-legacy-file-format>

function vtkwrite( filename,dataType,varargin )
% VTKWRITE Writes 3D Matlab array into VTK file format.
%  vtkwrite(filename,'structured_grid',x,y,z,'vectors',title,u,v,w) writes
%  a structured 3D vector data into VTK file, with name specified by the string
%  filename. (u,v,w) are the vector components at the points (x,y,z). x,y,z
%  should be 3-D matrices like those generated by meshgrid, where
%  point(ijk) is specified by x(i,j,k), y(i,j,k) and z(i,j,k).
%  The matrices x,y,z,u,v,w must all be the same size and contain
%  corrresponding position and vector component. The string title specifies
%  the name of the vector field to be saved. 
%
%  vtkwrite(filename,'structured_grid',x,y,z,'scalars',title,r) writes a 3D
%  scalar data into VTK file whose name is specified by the string
%  filename. r is the scalar value at the points (x,y,z). The matrices
%  x,y,z,r must all be the same size and contain the corresponding position
%  and scalar values. 
%
%  vtkwrite(filename,'structured_grid',x,y,z,'vectors',title1,u,v,w,'scalars',
%  title2,r) writes a 3D structured grid that contains both vector and scalar values.
%  x,y,z,u,v,w,r must all be the same size and contain the corresponding
%  positon, vector and scalar values.
%
%  vtkwrite(filename,'unstructured_grid',x,y,z,'vectors',title1,u,v,w,'scalars',
%  title2,r) writes a 3D unstructured grid that contains both vector and scalar values.
%  x,y,z,u,v,w,r must all be the same size and contain the corresponding
%  positon, vector and scalar values.
%  
%  vtkwrite(filename,'polydata','lines',x,y,z) exports a 3D line where
%  x,y,z are coordinates of the points that make the line. x,y,z are
%  vectors containing the coordinates of points of the line, where point(n)
%  is specified by x(n),y(n) and z(n).
%
%  vtkwrite(filename,'polydata','lines',x,y,z,'Precision',n) allows you to
%  specify precision of the exported number up to n digits after decimal
%  point. Default precision is 3 digits. 
%
%  vtkwrite(filename,'polydata','triangle',x,y,z,tri) exports a list of
%  triangles where x,y,z are the coordinates of the points and tri is an
%  m*3 matrix whose rows denote the points of the individual triangles.
%
%  vtkwrite(filename,'polydata','tetrahedron',x,y,z,tetra) exports a list
%  of tetrahedrons where x,y,z are the coordinates of the points
%  and tetra is an m*4 matrix whose rows denote the points of individual
%  tetrahedrons. 
%  
%  vtkwrite('execute','polydata','lines',x,y,z) will save data with default
%  filename ''matlab_export.vtk' and automatically loads data into
%  ParaView. 
%  
%  Version 2.2.1
%  Chaoyuan Yeh, 2015
%  Codes are modified from William Thielicke and David Gingras's submission.    

if strcmpi(filename,'execute'), filename = 'matlab_export.vtk'; end
fid = fopen(filename, 'w'); 
% VTK files contain five major parts
% 1. VTK DataFile Version
fprintf(fid, '# vtk DataFile Version 2.0\n');
% 2. Title
fprintf(fid, 'VTK from Matlab\n');
switch upper(dataType)
    case {'STRUCTURED_GRID','UNSTRUCTURED_GRID'}
        % 3. The format data proper is saved in (ASCII or Binary). Use
        % fprintf to write data in the case of ASCII and fwrite for binary.
        if numel(varargin)<6, error('Not enough input arguments'); end
        fprintf(fid, 'BINARY\n');
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
        if sum(size(x)==size(y) & size(y)==size(z))~=length(size(x))
            error('Input dimesions do not match')
        end
        n_elements = numel(x);
        % 4. Type of Dataset ( can be STRUCTURED_POINTS, STRUCTURED_GRID,
        % UNSTRUCTURED_GRID, POLYDATA, RECTILINEAR_GRID or FIELD )
        % This part, dataset structure, begins with a line containing the
        % keyword 'DATASET' followed by a keyword describing the type of dataset.
        % Then the geomettry part describes geometry and topology of the dataset.
        if strcmpi(dataType,'STRUCTURED_GRID')
            fprintf(fid, 'DATASET STRUCTURED_GRID\n');
            fprintf(fid, ['DIMENSIONS ' num2str(size(x,2)) ' ' num2str(size(x,1)) ' ' num2str(size(x,3)) '\n']);
        else
            fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
        end
        fprintf(fid, ['POINTS ' num2str(n_elements) ' float\n']);
        fwrite(fid, [x(:)';y(:)';z(:)'],'float','b');
        % 5.This final part describe the dataset attributes and begins with the
        % keywords 'POINT_DATA' or 'CELL_DATA', followed by an integer number
        % specifying the number of points of cells. Other keyword/data combination
        % then define the actual dataset attribute values.
        fprintf(fid, ['\nPOINT_DATA ' num2str(n_elements)]);
        % Parse remaining argument.
        vidx = find(strcmpi(varargin,'VECTORS'));
        sidx = find(strcmpi(varargin,'SCALARS'));
        if vidx~=0
            for ii = 1:length(vidx)
                title = varargin{vidx(ii)+1};
                % Data enteries begin with a keyword specifying data type
                % and numeric format.
                fprintf(fid, ['\nVECTORS ', title,' float\n']);
                fwrite(fid, [reshape(varargin{vidx(ii)+2},1,n_elements);...
                reshape(varargin{vidx(ii)+3},1,n_elements);...
                reshape(varargin{vidx(ii)+4},1,n_elements)],'float','b');
            end
        end
        if sidx~=0
            for ii = 1:length(sidx)
                title = varargin{sidx(ii)+1};
                fprintf(fid, ['\nSCALARS ', title,' float\n']);
                fprintf(fid, 'LOOKUP_TABLE default\n');
                fwrite(fid, reshape(varargin{sidx(ii)+2},1,n_elements),'float','b');
            end
        end
        
    case 'POLYDATA'
        fprintf(fid, 'ASCII\n');
        if numel(varargin)<4, error('Not enough input arguments'); end
        x = varargin{2}(:);
        y = varargin{3}(:);
        z = varargin{4}(:);
        if numel(varargin)<4, error('Not enough input arguments'); end
        if sum(size(x)==size(y) & size(y)==size(z))~= length(size(x))
            error('Input dimesions do not match')
        end
        n_elements = numel(x);
        fprintf(fid, 'DATASET POLYDATA\n');
        if mod(n_elements,3)==1
            x(n_elements+1:n_elements+2,1)=[0;0];
            y(n_elements+1:n_elements+2,1)=[0;0];
            z(n_elements+1:n_elements+2,1)=[0;0];
        elseif mod(n_elements,3)==2
            x(n_elements+1,1)=0;
            y(n_elements+1,1)=0;
            z(n_elements+1,1)=0;
        end
        nbpoint = numel(x);
        fprintf(fid, ['POINTS ' num2str(nbpoint) ' float\n']);
        precision = '3';
        if any(strcmpi(varargin,'PRECISION'))
            precision = num2str(uint8(varargin{find(strcmpi(varargin,'PRECISION'))+1}));
            if str2double(precision) < 0, error('Invalid precision spec.');end
        end
        spec = [repmat(['%0.',precision,'f '],1,9),'\n'];
        fprintf(fid,spec, [x(1:3:end-2) y(1:3:end-2) z(1:3:end-2) ...
             x(2:3:end-1) y(2:3:end-1) z(2:3:end-1) ...
             x(3:3:end) y(3:3:end) z(3:3:end)]');
        switch upper(varargin{1})
            case 'LINES'
                if mod(n_elements,2)==0
                    nbLine = 2*n_elements-2;
                else
                    nbLine = 2*(n_elements-1);
                end
                conn1 = zeros(nbLine,1);
                conn2 = zeros(nbLine,1);
                conn2(1:nbLine/2) = 1:nbLine/2;
                conn1(1:nbLine/2) = conn2(1:nbLine/2)-1;
                conn1(nbLine/2+1:end) = 1:nbLine/2;
                conn2(nbLine/2+1:end) = conn1(nbLine/2+1:end)-1;
                fprintf(fid,'\nLINES %d %d\n',nbLine,3*nbLine);
                fprintf(fid,'2 %d %d\n',[conn1';conn2']);
            case 'TRIANGLE'
                ntri = length(varargin{5});
                fprintf(fid,'\nPOLYGONS %d %d\n',ntri,4*ntri);
                fprintf(fid,'3 %d %d %d\n',(varargin{5}-1)');
            case 'TETRAHEDRON'
                ntetra = length(varargin{5});
                fprintf(fid,'\nPOLYGONS %d %d\n',ntetra,5*ntetra);
                fprintf(fid,'4 %d %d %d %d\n',(varargin{5}-1)');
        end     
end
fclose(fid);
if strcmpi(filename,'matlab_export.vtk')
    switch computer
        case {'PCWIN','PCWIN64'}
            !paraview.exe --data='matlab_export.vtk' &
            % Exclamation point character is a shell escape, the rest of the
            % input line will be sent to operating system. It can not take
            % variables, though. The & at the end of line will return control to 
            % Matlab even when the outside process is still running. 
        case {'GLNXA64','MACI64'}
            !paraview --data='matlab_export.vtk' &
    end
end


