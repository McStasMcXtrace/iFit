function out = openpdb(filename)
%OPENPDB Open a Protein Data Bank file, display its structure factor
%        and set the 'ans' variable to an iData object with its content

out = load(iData,filename, 'PDB');
subplot(out);

if ~isdeployed
  assignin('base','ans',out);
  ans = out
end
