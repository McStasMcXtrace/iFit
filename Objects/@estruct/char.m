function d = char(self)

  % build the output string
  d = '';
  for index=1:numel(self)
    s = self(index);
    if length(self) > 1
      d = [ d sprintf('%5i ',index) ];                        % index
    end
    if isempty(s.Tag)
      d = [ d sprintf('%9s ','<nul>') ];                      % Tag
    else
      d = [ d sprintf('%9s ',s.Tag) ];
    end
    d = [ d sprintf('%11s ', [ mat2str(size(s)) ]) ];  % size
    
    t = cellstr(s.Name); t = strtrim(t{1}); t(~isstrprop(t,'print') | t=='\' | t=='%')='';
    if length(t) > 31, t = [ t(1:27) '...' ]; end             % object.Name
    
    ts = title(s); if isempty(ts), ts = getaxis(s, '0'); end
    if ~ischar(ts), ts = ''; end
    t = [ t ' "' ts '"' ]; t = strtrim(t); t(~isstrprop(t,'print') | t=='\')='';
    if length(t) > 41, t = [ t(1:37) '..."'  ]; end           % title(Signal)
    d = [ d sprintf('%43s ', [ '''' t '''' ]) ];
    
    h = cellstr(s.Command); h = strtrim(h{end}); h(~isstrprop(h,'print') | h=='\')='';
    if length(h) > 23, h = [ h(1:20) '...' ]; end             % last command
    d = [ d sprintf('%s ', h) ];
    
    if ~isempty(s.Label) && ~isempty(s.DisplayName)
      h = cellstr(s.Label); h = strtrim(h{1}); h(~isstrprop(h,'print') | h=='\')='';
      if length(h) > 13, h = [ h(1:10) ]; end                 % Label/DisplayName
      d = [ d sprintf('%s', h) ];
      h = cellstr(s.DisplayName); h = strtrim(h{1}); h(~isstrprop(h,'print') | h=='\')='';
      if length(h) > 13, h = [ h(1:10) '...' ]; end           %
      d = [ d sprintf('/%s', h) ];
    elseif ~isempty(s.Label)
      h = cellstr(s.Label); h = strtrim(h{1}); h(~isstrprop(h,'print') | h=='\')='';
      if length(h) > 18, h = [ h(1:15) '...' ]; end           % Label
      d = [ d sprintf('%s', h) ];
    elseif ~isempty(s.DisplayName)
      h = cellstr(s.DisplayName); h = strtrim(h{1}); h(~isstrprop(h,'print') | h=='\')='';
      if length(h) > 18, h = [ h(1:15) '...' ]; end           % DisplayName
      d = [ d sprintf('%s', h) ];
    end

    d = [ d sprintf('\n') ];

  end
