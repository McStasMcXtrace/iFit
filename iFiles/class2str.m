function str=class2str(this, data)
% class2str(this,data) Create a string [ 'this = data;' ]
%   This function creates a string containing Matlab code describing a variable.
%
% input arguments:
%   this: string containg the name of the object to describe
%   data: any data set (struct, array, cell, iData, char)
%
% output variables:
%   str: string which contains a function code to generate the data.
%
% example: str=class2str('this', struct('a',1,'b','a comment','c',{});
%          
% See also: mat2str, num2str, eval, sprintf
%
% Part of: iFiles utilities (ILL library)
% Author:  E. Farhi <farhi@ill.fr>. June, 2007.

if nargin == 1
  data = this;
  if isempty(inputname(1)), this = [ class(data) '_str' ];
  else this = inputname(1); end
end

NL = sprintf('\n');
if ischar(data)
  str = [ this ' = ''' class2str_validstr(data) ''';' NL ];
elseif isa(data, 'iData')
  str = [ '% ' this ' (' class(data) ') size ' num2str(size(data)) NL ];
  str = [ str class2str(this, struct(data)) ];
  str = [ str NL '% handling of iData objects -------------------------------------' NL 'if ~exist(''iData''), return; end' NL ];
  str = [ str this '_s=' this '; ' this ' = rmfield(' this ',''Alias''); ' this ' = iData(' this ');' NL ...
         'setalias(' this ', ' this '_s.Alias.Names, ' this '_s.Alias.Values, ' this '_s.Alias.Labels);' NL ... 
         'if ~isempty(' this '_s.Alias.Axis)' NL ...
         '  setaxis('  this ', mat2str(1:length(' this '_s.Alias.Axis)), ' this '_s.Alias.Axis);' NL ...
         'end' NL ...
         '% end of iData ' this NL ];
elseif isnumeric(data) | islogical(data)
  str = [ '% ' this ' numeric (' class(data) ') size ' num2str(size(data)) NL ...
          this ' = ' mat2str(data(:)) ';' NL ];
  if prod(size(data)) > 1
    str = [ str this ' = reshape(' this ', [' num2str(size(data)) ']);' NL ];
  end
elseif isstruct(data)
  f = fieldnames(data);
  str = [ '% ' this ' (' class(data) ') length ' num2str(length(f)) NL ];
  for index=1:length(f)
    str = [ str class2str([ this '.' f{index} ], getfield(data, f{index})) ];
  end
  str = [ str '% end of struct ' this NL ];
elseif iscellstr(data)
  str = [ '% ' this ' (' class(data) 'str) size ' mat2str(size(data)) NL ...
          this ' = { ...' NL ];
  for index=1:length(data(:))
    str = [ str '  ''' class2str_validstr(data{index}) '''' ];
    if index < length(data(:)), str = [ str ', ' ]; end
    str = [ str ' ...' NL ];
  end
  str = [ str '}; ' NL ];
  if prod(size(data)) > 1
    str = [ str this ' = reshape(' this ', [' mat2str(size(data)) ']);' NL ];
  end
  str = [ str '% end of cellstr ' this NL ];
elseif iscell(data)
  str = [ '% ' this class(data) ' size ' mat2str(size(data)) NL ...
          this ' = cell(' mat2str(size(data)) ');' NL ];
  for index=1:length(data(:))
    str = [ str class2str([ this '{' num2str(index) '}' ], data{index}) NL ];
  end
  if prod(size(data)) > 1
    str = [ str this ' = reshape(' this ', [' mat2str(size(data)) ']);' NL ];
  end
  str = [ str '% end of cell ' this NL ];
elseif isa(data, 'function_handle')
  iData_private_warning(mfilename,  'can not save function handles. Skipping.');
else
  try
    % other class
    str = [ '% ' this ' (' class(data) ') size ' num2str(size(data)) NL ];
    str = [ str class2str(this, struct(data)) ];
    str = [ str '% end of object ' this NL ];
  catch
    iData_private_warning(mfilename,[ 'can not save ' class(data) '. Skipping.' ]);
  end
end

function str=class2str_validstr(str)
index = find(str < 32 | str > 127);
str(index) = ' ';
str=strrep(str, '''', '''''');
