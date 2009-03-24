function s = commandhistory(a)
% commandhistory(s) : show the command history of iData object
%
%   @iData/commandhistory shows the list of commands that have been used to
%   produce the current iData object.
%
% input:  s: object or array (iData)
% output: b: command history (char/cell)
% ex:     b=commandhistory(a);
%
% Version: $Revision: 1.1 $
% See also iData, iData/disp, iData/display

% handle input iData arrays
if length(a(:)) > 1
  s = {};
  for index=1:length(a(:))
    s{index} = commandhistory(a(index));
  end
  s = reshape(s, size(a));
  return
end

s = a.Command;
