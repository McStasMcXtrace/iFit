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
% Version: $Revision: 1.4 $
% See also iData, iData/disp, iData/display

% handle input iData arrays
if numel(a) > 1
  s = cell(size(a));
  parfor index=1:numel(a)
    s{index} = commandhistory(a(index));
  end
  return
end

s = a.Command;
if nargout == 0
  T   = a.Title; if iscell(T), T=T{1}; end
  T   = regexprep(T,'\s+',' '); % remove duplicated spaces
  selection = listdlg('ListString', s, 'ListSize',[400 300],'Name', T ,'PromptString', char(a)); 
  if ~isempty(selection)
      s=s{selection};
  end
end
