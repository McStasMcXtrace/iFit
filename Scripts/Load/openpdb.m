function out = openpdb(filename)
%OPENPDB Open a Protein Data Bank file, display its structure factor
%        and set the 'ans' variable to an iData object with its content
% (c) E.Farhi, ILL. License: EUPL.

if ~isa(filename,'iData')
  out = iData(filename,'PDB');  % with post-processing
else
  out = filename;
end
clear filename;

if numel(out) > 1
  % handle input iData arrays
  for index=1:numel(out)
    out(index) = feval(mfilename, out(index));
  end
end

if ~nargout
  figure; subplot(out);
  
  if ~isdeployed
    assignin('base','ans',out);
    ans = out
  end
end

