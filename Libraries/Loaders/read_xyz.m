% READ_XYZ (read a Molecule XYZ file)
%
%   Usage: s = read_xyz(filename)
%                    
% References:
%   https://en.wikipedia.org/wiki/XYZ_file_format
%   http://www.mathworks.com/matlabcentral/fileexchange/18233-nanovis--molecular-visualizer
%
% See also: read_poscar, read_pdb, read_mrc

function s = read_xyz(filename)

  s = [];
  if nargin == 0, return; end
  
  fid = fopen(filename,'r');


  if fid==-1
     return
  end

  % from http://www.mathworks.com/matlabcentral/fileexchange/18233-nanovis--molecular-visualizer/content/NanoVIs/NanoVis.m
  %Read a "xyz" file
    
  natom = fscanf(fid,'%i\n');
 
  fscanf(fid,'%i\n');
  data=zeros(natom,4);
  properties=[];
  species={};
 
  for k1=1:natom
     
    sp = fscanf(fid,'%c',2);
    species{end+1}=sp;       
      
    %Read coordinates
    data(k1,1)=fscanf(fid,'%f',1);  % XYZ
    data(k1,2)=fscanf(fid,'%f',1);
    data(k1,3)=fscanf(fid,'%f\n',1); 
    
    %read valence
    % valence(k1)=round(properties(data(k1,1),3));
     
  end 
  
  fclose(fid);
  
  % assemble the 'geometry' structure in the style returned by vasplab/import_poscar
  % http://www.mathworks.com/matlabcentral/fileexchange/36836-vasplab/content/vasplab/import_poscar.m
  s.filename = filename;
  s.comment  = comment;
  s.symbols  = species;
  s.coords   = data;
  s.selective = 0;
  s.atomcount = natom;
  
