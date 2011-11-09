function out = opensqw(filename)
%OPENSQW Open a McStas Sqw file (isotropic dynamic stucture factor, text file) 
%        display it and set the 'ans' variable to an iData object with its content

out = iData(filename);
plot(out);

if ~isdeployed
  assignin('base','ans',out);
end
