function out = openyaml(filename)
%OPENYAML Open a YAML/JSON File, display it
%        and set the 'ans' variable to an iData object with its content

out = load(iData,filename, 'YAML');
subplot(out);

if ~isdeployed
  assignin('base','ans',out);
  ans = out
end
