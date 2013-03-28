function out = openstl(filename)
%OPENSTL Open an STL/SLP 3D ascii stereolithography data file, display it
%        and set the 'ans' variable to an iData object with its content

out = load(iData,filename, 'STL');
subplot(out);

if ~isdeployed
  assignin('base','ans',out);
  ans = out
end
