function a=load_mcstas_sim(a0)
% function a=load_mcstas_sim(a0)
%
% Returns an iData style dataset from a McStas sim file
%
% Version: $Revision: 1.7 $
% See also: iData/load, iLoad, save, iData/saveas

% handle input iData arrays
if length(a0(:)) > 1
  a = [];
  for index=1:length(a0(:))
    a = [ a feval(mfilename, a0(index)) ];
  end
  return
end

if isempty(findstr(a0,'McStas'))
  warning([ mfilename ': The loaded data set ' a0.Tag ' from ' a0.Source ' is not a McStas data format.' ]);
  a = a0;
  return
end

% Find filename fields in sim struct:
filenames = findstr(a0,'filename');
dirname   = fileparts(a0.Source);
a=[];
if length(filenames(:)) > 0
  % This is a McStas 'overview' plot
  for j=1:length(filenames(:))
    filename = filenames{j};
    filename = strrep(filename,';','');
    filename(1:length('filename: '))='';
    filename(findstr(' ',filename):length(filename))='';
    if isempty(filename), continue; end
    filename = fullfile(dirname,filename);
    b = iData(filename);
    a = [a b];
  end
else
  % This is a .sim from a scan
  filename = 'mcstas.dat';
  filename = fullfile(dirname,filename);
  a = iData(filename);
end
