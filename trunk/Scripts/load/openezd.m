function out = openezd(filename)
%OPENEZD Open an EZD electronic density map File, display it
%        and set the 'ans' variable to an iData object with its content

out = load(iData,filename, 'PLY');
subplot(out);

if ~isdeployed
  assignin('base','ans',out);
  ans = out
end
