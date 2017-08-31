function [v,m] = mccode_run(m)
% run a McCode instrument model and request instrument parameters in a dialogue
% 
% You can change the simulation options and parameters from the dialogue shown.
% Then, close the window or select the 'OK' context menu item (right-click) to
% start the McDode calculation.
% The numeric input instrument parameters can be given as vectors (between [ ... ])
% in order to  start a scan calculation.
%
% To abort, select the CANCEL context menu item (right-click), or press Ctrl-C
% during the calculation.

% the 4 first items are simulation options
for f={'ncount','gravitation','mpi','seed'}; options.(f{1}) = m.UserData.options.(f{1}); end
% then we have numel(m.Parameters) scalar parameters
p1 = cell2struct(num2cell(m.ParameterValues(:)), strtok(m.Parameters(:)), 1);
% and finally numel(fieldnames(m.UserData.Parameters_Constant)) others
p2 = m.UserData.Parameters_Constant;
dlg.Name = m.Name;
dlg.TooltipString = m.Name;
dlg.ColumnFormat = 'char';  % this way we can change parameter values from scalar to vectors
answer = inputdlg(cat(options, p1, p2), dlg);

% check if we have cancelled the run action
if isempty(answer), v=[]; return; end

% transfer the answer into the model properties
for f={'ncount','gravitation','mpi','seed'}; m.UserData.options.(f{1}) = answer.(f{1}); end
for f=m.Parameters(:)'; 
  m.(f{1}) = answer.(f{1}); 
  p.(f{1}) = answer.(f{1}); % can also be a vector (scan)
end
if ~isempty(m.UserData.Parameters_Constant)
  for f=fieldnames(m.UserData.Parameters_Constant)
    m.UserData.Parameters_Constant.(f{1}) = answer.(f{1});
  end
end

if ~isempty(inputname(1)) && ~isempty(m)
  assignin('caller', inputname(1), m);
end

% run it and get back result
v = iData(m, p, nan); % we want the raw data



