function str=class2str(this, data, flat)
% class2str(this,data) Create a string [ 'this = data;' ]
%   This function creates a string containing Matlab code describing a variable.
%   The syntax class2str(this, data, flat) creates a flat text with commented 
%   data blocks, which is not an m-file, but rather a Linux-style config file.
%
% input arguments:
%   this: string containg the name of the object to describe
%   data: any data set (struct, array, cell, iData, char)
%   flat: optinal argument to generate a flat file
%
% output variables:
%   str: string which contains a function code to generate the data.
%
% example: str=class2str('this', struct('a',1,'b','a comment','c',{});
%          
% See also: mat2str, num2str, eval, sprintf
%
% Part of: iFiles utilities (ILL library)
% Author:  E. Farhi <farhi@ill.fr>. $Revision: 1.6 $

if nargin == 1
  data = this;
  if isempty(inputname(1)), this = [ class(data) '_str' ];
  else this = inputname(1); end
end

if nargin == 3,
  str = class2str_flat(this, data);
else
  str = class2str_m(this, data);
end

return

% ------------------------------------------------------------------------------
function str = class2str_m(this, data)
% function to create an m-file string

NL = sprintf('\n');
if ischar(data)
  str = [ this ' = ''' class2str_validstr(data) ''';' NL ];
elseif (isa(data, 'iData') | isstruct(data)) & length(data) > 1
  str = [ '% ' this ' (' class(data) ') array size ' mat2str(size(data)) NL ];
  for index=1:length(data(:))
    str = [ str class2str([ this '(' num2str(index) ')' ], data(index)) NL ];
  end
  str = [ str this ' = reshape(' this ', [' num2str(size(data)) ']);' NL ...
          '% end of ' class(data) ' array ' this NL ];
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
  str = [ '% ' this ' function (' class(data) ')' NL ...
          this ' = ' func2str(data(:)) ';' NL ];
else
  try
    % other class
    str = [ '% ' this ' (' class(data) ') size ' num2str(size(data)) NL ];
    str = [ str class2str(this, struct(data)) ];
    str = [ str '% end of object ' this NL ];
  catch
    warning([ mfilename ': can not save ' this ' ' class(data) '. Skipping.' ]);
  end
end

% ------------------------------------------------------------------------------
function str = class2str_flat(this, data)
% function to create a flat file string

str = '';
NL  = sprintf('\n');
if isempty(data), return; end
if ischar(data)
  str = [ '# ' this ': ' class2str_validstr(data) NL ];
elseif (isa(data, 'iData') | isstruct(data)) & length(data) > 1
  for index=1:length(data(:))
    str = [ str class2str([ this '(' num2str(index) ')' ], data(index), 'flat') NL ];
  end
elseif isstruct(data)
  f = fieldnames(data);
  %str = [ '# ' class(data) ' length ' num2str(length(f)) ': ' this NL ];
  str = '';
  for index=1:length(f)
    str = [ str class2str([ this ' ' f{index} ], getfield(data, f{index}), 'flat') ];
  end
elseif iscellstr(data)
  % str = [ '# ' class(data) 'str size ' mat2str(size(data)) ': ' this NL ];
  str = '';
  for index=1:length(data(:))
    str = [ str class2str(this,data{index},'flat') ];
  end
elseif iscell(data)
  % str = [ '# ' class(data) ' size ' mat2str(size(data)) ': ' this NL ];
  str = '';
  for index=1:length(data(:))
    str = [ str class2str([ this '{' num2str(index) '}' ], data{index}, 'flat') NL ];
  end
elseif isa(data, 'function_handle')
  str = [ '# ' this ' function: ' func2str(data(:)) NL ];
elseif (isnumeric(data) & ~isa(data,'iData')) | islogical(data)
  str = num2str(data);
  str(:,end+1) = sprintf('\n');
  str = str';
  str = str(:)';
  str = [ '# ' class(data) ' size ' num2str(size(data)) ': ' this NL ...
          str ];
else
  try
    % other class
    str = [ '# ' class(data) ' size ' num2str(size(data)) ': ' this NL ];
    str = [ str class2str(this, struct(data),'flat') ];
  catch
    warning([ mfilename ': can not save ' this ' ' class(data) '. Skipping.' ]);
  end
end

% ------------------------------------------------------------------------------

function str=class2str_validstr(str)
% validate a string as a single line
str=strrep(str(:)', sprintf('\n'), ';');
index = find(str < 32 | str > 127);
str(index) = ' ';
str=strrep(str, '''', '''''');
