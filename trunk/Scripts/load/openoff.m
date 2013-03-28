function out = openoff(filename)
%OPENOFF Open an OFF 3D ascii File, display it
%        and set the 'ans' variable to an iData object with its content

out = load(iData,filename, 'OFF');
subplot(out);

if ~isdeployed
  assignin('base','ans',out);
  ans = out
end
