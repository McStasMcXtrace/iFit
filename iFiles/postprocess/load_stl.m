function b=load_stl(a)
% function a=load_stl(a)
%
% Returns an iData style dataset from a STL file (ascii or binary)
%
% Version: $Revision: 1.1 $
% See also: iData/load, iLoad, save, iData/saveas

% handle input iData arrays
if length(a(:)) > 1
  for index=1:length(a(:))
    a(index) = feval(mfilename, a(index));
  end
  return
end

if ~isfield(a.Data, 'vertices')
  warning([ mfilename ': The loaded data set ' a.Tag ' from ' a.Source ' is not a STL data format.' ]); 
  return
end

a.Data.Signal=ones(size(a.Data.vertices, 1),1);
setalias(a, 'Signal', 'Data.Signal');
setalias(a, 'X', 'Data.vertices(:,1)');
setalias(a, 'Y', 'Data.vertices(:,2)');
setalias(a, 'Z', 'Data.vertices(:,3)');

setaxis(a, 1, 'X'); setaxis(a, 2, 'Y'); setaxis(a, 3, 'Z');


b=a;
