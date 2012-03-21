function b=load_stl(a)
% function a=load_stl(a)
%
% Returns an iData style dataset from a STL file (ascii or binary)
%
% Version: $Revision: 1.2 $
% See also: iData/load, iLoad, save, iData/saveas

% handle input iData arrays
if length(a(:)) > 1
  b = [];
  for index=1:length(a(:))
    b = [ b feval(mfilename, a(index)) ];
  end
  return
end

if isfield(a.Data.MetaData, 'OFF')
  % this is an OFF format file read by looktxt
  nvf=a.Data.MetaData.OFF;      % 'NVertices  NFaces  NEdges'
  nv=nvf(1,1); nf=nvf(1,2); % indices start at 0 in OFF
  a.Data.face = nf; a.Data.vertex = nv;
  a.Data.MetaData.OFF = [nv nf];
  nvf=nvf(2:end,:);
elseif strncmpi(a.Format, 'ply',3)
  nf = a.Data.face;
  nv = a.Data.vertex;
  % the data block is Data.end_header in all cases
  nvf= a.Data.end_header;
  a.Data.end_header = [];
end

if isfield(a.Data.MetaData, 'OFF') || strncmpi(a.Format, 'ply',3)
  if size(nvf,1) <= nv+nf  % only contains vertices: an other block gives the faces...
    a.Data.vertices=nvf(1:nv,:);
    [match,type,n] = findfield(a);
    % sort fields by size
    [n,sorti]=sort(n,'descend');
    for index=1:length(n)
      if ~strcmp(type,'double'), continue; end  % only get block with numeric data
      if n(index) == numel(nvf), continue; end   % this block has already been found (vertices)
      f=get(a,match{sorti(index)});
      if isfield(a.Data.MetaData, 'OFF'), f=f+1; end  % indices start at 0 in OFF
      a.Data.faces = f; clear f;  
      set(a,match{sorti(index)},[]);
      break
    end
  else
    a.Data.vertices=nvf(1:nv,:);
    a.Data.faces   =nvf(nv:end,:);
  end
end  
if ~isfield(a.Data, 'vertices')
  warning([ mfilename ': The loaded data set ' a.Tag ' from ' a.Source ' is not a STL/SLP/OFF/PLY data format.' ]); 
  b = a;
  return
end

a.Data.Signal=ones(size(a.Data.vertices, 1),1);
setalias(a, 'Signal', 'Data.Signal');
setalias(a, 'X', 'Data.vertices(:,1)');
setalias(a, 'Y', 'Data.vertices(:,2)');
setalias(a, 'Z', 'Data.vertices(:,3)');

setaxis(a, 1, 'X'); setaxis(a, 2, 'Y'); setaxis(a, 3, 'Z');


b=a;
