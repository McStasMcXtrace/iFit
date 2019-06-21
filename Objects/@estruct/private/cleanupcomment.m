function str = cleanupcomment(comment,option)
% CLEANUPCOMMENT clean a string from EOL and duplicate spaces.
%   S = CLEANUPCOMMENT(C) removes non printable characters and duplicate spaces.
%
%   S = CLEANUPCOMMENT(C,'long') same as above, but keeps duplicate spaces. 

  if nargin < 2, option='short'; end
  
  % Replace linefeeds and carriage returns.
  str = strrep(comment, '\', '/');
  isprint = isstrprop(str, 'print');
  str(~isprint) = ' ';

  if any(strcmp(option,{'short','compact'}))
    % Replace all double spaces with single spaces.
    while (strfind(str, '  '))
	    str = strrep(str, '  ', ' ');
    end
  end

  % Remove leading and trailing space.
  str = strtrim(str);

