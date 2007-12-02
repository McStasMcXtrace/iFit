function a=mcsimload(a0)
% function a=mcsimload(a0)
%
% Returns an iData style dataset from a McStas sim file
%

% Find filename fields in sim struct:
filenames = findstr(a0,'filename');
dirname = fileparts(a0.Source);

a=[];

for j=1:length(filenames(:))
  filename = filenames{j};
  filename(1:length('filename: '))='';
  filename(findstr(' ',filename):length(filename))='';
  filename = fullfile(dirname,filename);
  b = iData(filename);
  a = [a b];
end