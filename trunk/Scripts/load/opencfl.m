function out = opencfl(filename)
%OPENCFL Open a CFL FullProf crystallography file, display it
%        and set the 'ans' variable to an iData object with its content
%        (Required cif2hkl to have been compiled with gfortran, see Install page).

out = load(iData,filename, 'CFL');
subplot(out);

if ~isdeployed
  assignin('base','ans',out);
  ans = out
end
