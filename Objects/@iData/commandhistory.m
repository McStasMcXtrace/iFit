function [s,fig] = commandhistory(a, fig)
% commandhistory(s) : show the command history of iData object
%
%   @iData/commandhistory shows the list of commands that have been used to
%   produce the current iData object.
%
% input:  s: object or array (iData)
% output: b: command history (char/cell)
% ex:     b=commandhistory(a);
%
% Version: $Date$
% See also iData, iData/disp, iData/display

if nargin == 2 && ishandle(fig)
  s=commandhistory_export(fig);
  return
end

% handle input iData arrays
if numel(a) > 1
  s = cell(size(a)); fig=s;
  parfor index=1:numel(a)
    [s{index},fig{index}] = commandhistory(a(index));
  end
  return
end

s = a.Command;
if nargout == 0 || nargout == 2
  T   = a.Title; if iscell(T), T=T{1}; end
  T   = regexprep(T,'\s+',' '); % remove duplicated spaces
  [fig] = listdlg_nonmodal('ListString', s, 'ListSize',[400 300], ...
    'Name', T , ...
    'PromptString', char(a), 'OKString','Save all to file...','CancelString','OK'); 

  ad = getappdata(fig,'ListDialogAppData__');
  ad.object = a;
  setappdata(fig,'ListDialogAppData__', ad);
    
  set(fig, 'closerequestfcn','commandhistory(iData, gcbf); delete(gcbf)');
  selection=[]; ok=0;
end
  
function s=commandhistory_export(fig)

  ad = getappdata(fig,'ListDialogAppData__');
  if strcmpi(ad.button, 'cancel'), s=[]; return; end
  s = get(ad.listbox,'string');
  if isempty(s), return; end
  a        = ad.object;
  
  % save all commands to a script file
  [filename, pathname] = uiputfile('*.m', 'Save commands to file', [ 'iFit_' a.Tag '_history' ]);
  if filename == 0, return; end % user press Cancel
  filename = fullfile(pathname, filename);
  [pathname, name, ext] = fileparts(filename); % get name without extension
  [fid, message] = fopen(filename, 'w+');
  if fid == -1    % invalid filename
    iData_private_warning(mfilename,[ 'Error opening file ' filename ' to save object command history' ]);
    disp(message)
    return
  end
  
  % assemble the header
  NL = sprintf('\n');
  str = [ '% Matlab script generated by iFit/iData/commandhistory' NL ...
          '% File: ' filename NL ...
          '%   To use/import data, type "' name '" at the matlab prompt.' NL ...
          '%   You will obtain an iData object (if you have iFit/iData installed) or a structure.' NL ...
          '%   This script may depend on other objects, which you have to create with similar scripts.' NL ...
          '% Original data: ' NL ...
          '%   class:    ' class(a) NL ...
          '%   variable: ' inputname(1) NL ...
          '%   tag:      ' a.Tag NL ...
          '%   label:    ' a.Label NL ...
          '%   source:   ' a.Source NL ... 
          '%' NL ...
          '% Matlab ' version ' m-file ' filename ' saved on ' datestr(now) NL ...
          '%' NL ];

  fprintf(fid, '%s', str);
  for index=1:length(s)
    fprintf(fid, '%s\n', s{index});
  end
  fprintf(fid, '%% End of script %s\n', filename);
  fclose(fid);

