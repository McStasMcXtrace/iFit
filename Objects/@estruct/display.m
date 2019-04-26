function d = display(s_in, name)
% d = display(s) : display iData object (from command line)
%
% @iData/display function to display iData object.
%   Used when no ';' sign follows a iData object in matlab.
% The return value may be catched as a string to display.
%
% input:  s: object or array (iData)
% output: d: string to display (char)
% ex:     'display(iData)' or 'iData'
%
% Version: $Date$ $Version$ $Author$
% See also iData, iData/disp, iData/get

% EF 27/07/00 creation
% EF 23/09/07 iData implementation

if nargin == 2 && ~isempty(name)
  iname = name;
elseif ~isempty(inputname(1))
  iname = inputname(1);
else
  iname = 'ans';
end

% build the header -----------------------------------------------------
d = [ sprintf('%s = ',iname) ];

if numel(s_in) > 1
  d = [ d sprintf(' array [%s]',num2str(size(s_in))) ];
end
if isdeployed || ~usejava('jvm') || ~usejava('desktop'), id=class(s_in);
else           id=[ '<a href="matlab:doc ' class(s_in) '">' class(s_in) '</a> ' ...
                    '(<a href="matlab:methods ' class(s_in) '">methods</a>,' ...
                    '<a href="matlab:doc(''' class(s_in) ''')">doc</a>,' ...
                    '<a href="matlab:figure;subplot(' iname ');">plot</a>,' ...
                    '<a href="matlab:figure;subplot(log(' iname '));">plot(log)</a>,' ...
                    '<a href="matlab:double(' iname ')">values</a>,' ...
                    '<a href="matlab:disp(' iname ');">more...</a>)' ];
end
if isvector(s_in) > 1, id = [ id ' list/event']; end
if length(s_in) == 0
  d = [ d sprintf(' %s object: empty\n',id) ];
else
  if numel(s_in) == 1
    d = [ d sprintf(' %s %iD object:\n\n', id, ndims(s_in)) ];
  else
    d = [ d sprintf(' %s object:\n\n', id) ];
  end
  if numel(s_in) > 1
    d = [ d sprintf('Index') ];
  end
  d = [ d sprintf('    [Tag] [Dimension]                                     [Title] [Last command]') ];
  if numel(s_in) > 1
    if any(~cellfun('isempty', get(s_in,'Label'))) || any(~cellfun('isempty', get(s_in,'DisplayName')))
      d = [ d '          [Label/DisplayName]' ];
    end
  else
    if ~isempty(get(s_in,'Label')) || ~isempty( get(s_in,'DisplayName'))
      d = [ d '          [Label/DisplayName]' ];
    end
  end
  d = [ d sprintf('\n') ];


  % add the objects representation
  d = [ d char(s_in) sprintf('\n') ];
end

if nargout == 0
  fprintf(1,d);
end