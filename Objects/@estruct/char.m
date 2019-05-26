function d = char(self, option)
% CHAR Convert object to character string.
%   S = CHAR(X) converts the object X into a character representation,
%   showing its Tag, Title/Name, last command, and Label.
%
%   S = CHAR(X, 'short') and S = CHAR(X, 'compact') produce a compact representation
%
%   Use CELLSTR to obtain an exact string representation which evaluation
%   rebuilds the object.
%
% Example: s=estruct(1:10); ischar(char(s))
% Version: $Date$ $Version$ $Author$
% See also: estruct, cellstr, double.

  if nargin < 2, option=''; end
  % build the output string
  d = '';
  for index=1:numel(self)
    s = self(index);
    if length(self) > 1, d1 = index;   else d1=[]; end        % index
    if isempty(s.Tag),   d2 = '<nul>'; else d2=s.Tag; end     % Tag
    d3 = mat2str(size(s));                                    % size

    t = cellstr(s.Name); t = strtrim(t{1}); t(~isstrprop(t,'print') | t=='\' | t=='%')='';
    if length(t) > 31, t = [ t(1:27) '...' ]; end             % Name

    ts = title(s);
    if ~ischar(ts), ts = ''; end
    t = [ t ' "' ts '"' ]; t = strtrim(t); t(~isstrprop(t,'print') | t=='\')='';
    if length(t) > 41, t = [ t(1:37) '..."'  ]; end           % title(Signal)
    d4 = [ '''' t '''' ];

    h = cellstr(s.Command); if ~isempty(h)
      h = strtrim(h{end}); h(~isstrprop(h,'print') | h=='\')='';
      if length(h) > 23, h = [ h(1:20) '...' ]; end             % last command
    end
    d5 = char(h);

    if ~isempty(s.Label)
      h = cellstr(s.Label); h = strtrim(h{1}); h(~isstrprop(h,'print') | h=='\')='';
      if length(h) > 18, h = [ h(1:15) '...' ]; end                 % Label/DisplayName
      d6 = h;
    else d6=[]; end
    if ~isempty(s.DisplayName)
      h = cellstr(s.DisplayName); h = strtrim(h{1}); h(~isstrprop(h,'print') | h=='\')='';
      if length(h) > 18, h = [ h(1:15) '...' ]; end           % DisplayName
      d7 = h;
    else d7 = []; end
    if ~isempty(d6) && ~isempty(d7)
      d6=[ d6 '/' ];
    end

    % build the final string
    if nargin == 1 || ~any(strcmp(option,{'short','compact'}))
      d = [ d cleanupcomment(sprintf('%5i %8s %11s %43s %23s %s%s', ...
        d1, d2, d3, d4, d5, d6, d7),'long') ];
    else % compact form
      d = [ d cleanupcomment(sprintf('%i %s %s %s %s %s%s', ...
        d1, d2, d3, d4, d5, d6, d7)) ];
    end
    if numel(self) > 1
      d = [ d sprintf('\n') ];
    end

  end
