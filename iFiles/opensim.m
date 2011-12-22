function out = opensim(filename)
%OPENSIM Open a McStas SIM file, display it
%        and set the 'ans' variable to an iData object with its content

out = iData(filename);
plot(out);

if ~isdeployed
  assignin('base','ans',out);
  ans = out
end
