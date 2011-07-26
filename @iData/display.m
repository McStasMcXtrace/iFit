function d = display(s_in)
% d = display(s) : display iData object (from command line)
%
% @iData/display function to display iData object.
%   Used when no ';' sign folows a iData object in matlab.
% The return value may be catched as a string to display.  
%
% input:  s: object or array (iData) 
% output: d: string to display (char)
% ex:     'display(iData)' or 'iData'
%
% Version: $Revision: 1.6 $
% See also iData, iData/disp, iData/get

% EF 27/07/00 creation
% EF 23/09/07 iData implementation

if ~isempty(inputname(1))
d = [ sprintf('%s = ',inputname(1)) ];
else
d = [ sprintf('%s = ','ans') ];
end
if length(s_in) > 1
  d = [ d sprintf(' array [%s]',num2str(size(s_in))) ];
end
if length(s_in) == 0
    d = [ d sprintf(' iData object: empty\n') ];
else
    d = [ d sprintf(' iData object:\n\n') ];
    s_in = s_in(:);
    if length(s_in) > 1
      d = [ d sprintf('Index ') ];
    end
    d = [ d sprintf('    [Tag] [Dimension]                                     [Title] [Last command]          [Label]\n') ];

    % now build the output string
    for index=1:length(s_in)
      s = s_in(index);
      if length(s_in) > 1
        d = [ d sprintf('%5i ',index) ];                        % index
      end
      if isempty(s.Tag)
        d = [ d sprintf('%9s ','<nul>') ];                      % Tag
      else
        d = [ d sprintf('%9s ',s.Tag) ];
      end
      d = [ d sprintf('%11s ', ['[' num2str(size(s)) ']' ]) ];  % size
      t = cellstr(s.Title); t = strtrim(t{1});
      if length(t) > 31, t = [ t(1:27) '...' ]; end             % object.title
      t = [ t ' "' title(s) '"' ];
      if length(t) > 41, t = [ t(1:37) '..."'  ]; end              % title(Signal)
      d = [ d sprintf('%43s ', [ '''' t '''' ]) ];
      h = cellstr(s.Command); h=deblank(h{end});
      if length(h) > 23, h = [ h(1:20) '...' ]; end             % last command
      d = [ d sprintf('%s ', h) ];
      h = cellstr(s.Label); h=deblank(h{end});
      if length(h) > 13, h = [ h(1:10) '...' ]; end             % Label
      d = [ d sprintf('%s ', h) ];

      d = [ d sprintf('\n') ];
    end
end



if nargout == 0
  fprintf(1,d);
end

