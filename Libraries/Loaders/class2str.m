function str=class2str(this, data, options)
% CLASS2STR This function creates a string containing Matlab code describing a variable.
%   CLASS2STR(this,data) Create an evaluable string [ 'this = data;' ]
%
%   CLASS2STR(data) is the same as above, and uses the variable name for
%   the label of the output.
%
%   CLASS2STR(..., 'no comments') same as above, but removes
%   comments from the output file.
%
%   CLASS2STR(..., 'eval') creates a compact evaluable string of
%   the initial data.
%
%   CLASS2STR(..., 'flat') creates a flat text with commented
%   data blocks, which is not an m-file, but rather a Linux-style config file.
%
%   CLASS2STR(..., 'short') strip-down the large elements in data
%   to reduce the total size of the representation. The output does not
%   contain the full initial object, but a reduced version.
%
%   Multiple keywords can be used, such as 'no comment, short'.
%
% Example: y=class2str(struct('a',1,'b','a string comment','c',{'cell'})); ischar(y)
%
%  data=struct('a',1,'b','a string comment','c',{'cell'});
%  str=class2str(data)          % produces a string/script that regenerates data
%  str=class2str('this', data)  % idem, but creates 'this' instead of 'data'.
%  str=class2str(data, 'flat')  % a kind of config file
%  str=class2str(data, 'eval')  % a compact evaluable string
%
% Version: $Date$ $Author$ $Author$
% See also: mat2str, num2str, eval, sprintf

  str=[];
  if     nargin == 0
    return
  elseif nargin == 1
    data    = this;
    this    = '';  % default root level name
    options = '';
  elseif nargin == 2
    if ~ischar(this) && ischar(data) % (data, options)
      options=data; data=this; this=inputname(1);
    else options = ''; end
  end
  if isempty(this), this = inputname(1); end
  if isempty(this), this = [ class(data) '_str' ]; end

  % handle a potential header
  nocomment = any(strfind(options, 'no comment'));
  NL = sprintf('\n');

  if any(strfind(options, 'flat'))
    if ~nocomment, str = [ '# created by iFit/class2str <ifit.mccode.org> on ' datestr(now) 'with option ' options NL ]; end
    str = [ str class2str_flat(this, data, options) ];
  elseif any(strfind(options, 'eval'))
    str = class2str_eval(data, options);
  else
    if ~nocomment, str = [ '% created by iFit/class2str <ifit.mccode.org> on ' datestr(now) ' with option ' options NL ]; end
    str = [ str class2str_m(this, data, options) ];
  end

% end class2str

% ------------------------------------------------------------------------------
function str = class2str_m(this, data, options)
% function to create an m-file string

  nocomment = any(strfind(options, 'no comment'));
  shorten   = any(strfind(options, 'short'));
  NL        = sprintf('\n');
  str       = '';


  if ischar(data)
    if shorten && numel(data) > 100, data = data'; data=[ transpose(data(1:97)) '...' ]; end
    str = [ this ' = ''' class2str_validstr(data) ''';' NL ];
  elseif (all(isobject(data)) || all(isstruct(data))) && numel(data) > 1
    if ~nocomment, str = [ '% ' this ' (' class(data) ') array size ' mat2str(size(data)) NL ]; end
    for index=1:numel(data)
      str = [ str class2str_m([ this '(' num2str(index) ')' ], data(index), options) ];
    end
    str = [ str this ' = reshape(' this ', [' num2str(size(data)) ']);' NL ];
    if ~nocomment, str = [ str  '% end of ' class(data) ' array ' this NL ]; end
  elseif isa(data, 'iData')
    if ~nocomment, str = [ '% ' this ' (' class(data) ') size ' num2str(size(data)) NL ]; end
    str = [ str class2str_m(this, struct(data), options) ];
    if ~nocomment, str = [ str NL '% handling of iData objects -------------------------------------' NL ]; end
    str = [ str 'if ~exist(''iData''), return; end' NL ];
    str = [ str this '_s=' this '; ' this ' = rmfield(' this ',''Alias''); ' this ' = iData(' this ');' NL ...
           'setalias(' this ', ' this '_s.Alias.Names, ' this '_s.Alias.Values, ' this '_s.Alias.Labels);' NL ...
           'if ~isempty(' this '_s.Alias.Axis)' NL ...
           '  setaxis('  this ', mat2str(1:length(' this '_s.Alias.Axis)), ' this '_s.Alias.Axis);' NL ...
           'end' NL ];
    if ~nocomment, str = [ str '% end of iData ' this NL ]; end
  elseif isobject(data)
    if ~nocomment, str = [ '% ' this ' (' class(data) ') size ' num2str(size(data)) NL ]; end
    str = [ str class2str_m(this, struct(data), options) ];
    if ~nocomment, str = [ str NL '% handling of objects -------------------------------------' NL ]; end
    str = [ str 'try; ' this ' = ' class(data) '(' this '); end' ];
    if ~nocomment, str = [ str '% end of object ' this NL ]; end
  elseif isnumeric(data) | islogical(data)
    if shorten && numel(data) > 100, data = data(1:100); end
    if ~nocomment, str = [ '% ' this ' numeric (' class(data) ') size ' num2str(size(data)) NL ]; end
    str = [ str this ' = ' mat2str(data(:)) ';' NL ];
    if numel(data) > 1
      str = [ str this ' = reshape(' this ', [' num2str(size(data)) ']);' NL ];
    end
  elseif isstruct(data) && ~isempty(data)
    f = fieldnames(data);
    if ~nocomment, str = [ '% ' this ' (' class(data) ') length ' num2str(length(f)) NL ]; end
    for index=1:length(f)
      if isempty(deblank(this))
        str = [ str class2str_m([ f{index} ], getfield(data, f{index}), options) ];
      else
        str = [ str class2str_m([ this '.' f{index} ], getfield(data, f{index}), options) ];
      end
    end
    if ~nocomment, str = [ str '% end of struct ' this NL ]; end
  elseif iscellstr(data)
    if ~nocomment, str = [ '% ' this ' (' class(data) 'str) size ' mat2str(size(data)) NL ]; end
    str = [ str this ' = { ...' NL ];
    for index=1:numel(data)
      str = [ str '  ''' class2str_validstr(data{index}) '''' ];
      if index < numel(data), str = [ str ', ' ]; end
      str = [ str ' ...' NL ];
    end
    str = [ str '}; ' NL ];
    if prod(size(data)) > 1
      str = [ str this ' = reshape(' this ', [' mat2str(size(data)) ']);' NL ];
    end
    if ~nocomment, str = [ str '% end of cellstr ' this NL ]; end
  elseif iscell(data)
    if ~nocomment, str = [ '% ' this class(data) ' size ' mat2str(size(data)) NL ]; end
    str = [ str this ' = cell(' mat2str(size(data)) ');' NL ];
    for index=1:numel(data)
      str = [ str class2str_m([ this '{' num2str(index) '}' ], data{index}, options) ];
    end
    if prod(size(data)) > 1
      str = [ str this ' = reshape(' this ', [' mat2str(size(data)) ']);' NL ];
    end
    if ~nocomment, str = [ str '% end of cell ' this NL ]; end
  elseif isa(data, 'function_handle')
    if ~nocomment, str = [ '% ' this ' function (' class(data) ')' NL ]; end
    str = [ str this ' = ' func2str(data(:)) ';' NL ];
  else
    try
      % other classes
      if ~nocomment, str = [ '% ' this ' (' class(data) ') size ' num2str(size(data)) NL ]; end
      str = [ str class2str_m(this, struct(data), options) ];
      str = [ str 'try; ' this ' = ' class(data) '(' this '); end' ];
      if ~nocomment, str = [ str '% end of ' class(data) ' ' this NL ]; end
    catch
      warning([ mfilename ': can not save ' this ' ' class(data) '. Skipping.' ]);
    end
  end
% end class2str_m

% ------------------------------------------------------------------------------
function str = class2str_flat(this, data, options)
% function to create a flat file string (no matlab code)

  nocomment = any(strfind(options, 'no comment'));
  shorten   = any(strfind(options, 'short'));
  NL        = sprintf('\n');
  str       = '';

  if isempty(data), return; end
  if ischar(data)
    if shorten && numel(data) > 100, data = data'; data=[ transpose(data(1:97)) '...' ]; end
    str = [ '# ' this ': ' class2str_validstr(data) NL ];
  elseif (isobject(data) || isstruct(data)) && length(data) > 1
    for index=1:numel(data)
      str = [ str class2str_flat([ this '(' num2str(index) ')' ], data(index), options) NL ];
    end
  elseif isstruct(data)
    f = fieldnames(data);
    %str = [ '# ' class(data) ' length ' num2str(length(f)) ': ' this NL ];
    str = '';
    for index=1:length(f)
      str = [ str class2str_flat([ this ' ' f{index} ], getfield(data, f{index}), options) ];
    end
  elseif iscellstr(data)
    % str = [ '# ' class(data) 'str size ' mat2str(size(data)) ': ' this NL ];
    str = '';
    for index=1:numel(data)
      str = [ str class2str_flat(this,data{index},options) ];
    end
  elseif iscell(data)
    % str = [ '# ' class(data) ' size ' mat2str(size(data)) ': ' this NL ];
    str = '';
    for index=1:numel(data)
      str = [ str class2str_flat([ this '{' num2str(index) '}' ], data{index}, options) NL ];
    end
  elseif isa(data, 'function_handle')
    str = [ '# ' this ' function: ' func2str(data(:)) NL ];
  elseif ~isobject(data) && (isnumeric(data) || islogical(data))
    if shorten && numel(data) > 100, data = data(1:100); end % starts with 1st column
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
      str = [ str class2str_flat(this, struct(data), options) ];
    catch
      warning([ mfilename ': can not save ' this ' ' class(data) '. Skipping.' ]);
    end
  end
% end class2str_flat

% ------------------------------------------------------------------------------
function str = class2str_eval(data, options)

  % can not handle comments as this is a single line
  shorten   = any(strfind(options, 'short'));
  NL        = sprintf('\n');
  str       = '';

  if ischar(data)
    if shorten && numel(data) > 100, data = data'; data=[ transpose(data(1:97)) '...' ]; end
    str = [ '''' class2str_validstr(data) '''' ];
  elseif numel(data) == 1 && ishandle(data) && ~isnumeric(data)
    % ignore, except for numeric (scalar)
  elseif isobject(data) && numel(data) == 1
    str = [ str class(data) '(' class2str_eval(struct(data), options) ')' ];
  elseif isnumeric(data)
    if shorten && numel(data) > 100, data = data(1:100); end % starts with 1st column
    if ndims(data) <= 2, str = mat2str(data); end
  elseif iscell(data)
    if isempty(data), str = [     '{ ' ]; end
    for index=1:numel(data)
      if index == 1, str = [     '{ ' class2str_eval(data{index}, options) ];
      else           str = [ str ', ' class2str_eval(data{index}, options) ]; end
    end
    str = [ str ' }' ];
    if size(data,1) > 1 && size(data,2) == 1, str = [ str '''' ]; end  % transpose
  elseif isa(data, 'function_handle') && numel(data) == 1
    str = func2str(data);
  elseif isstruct(data) && numel(data) == 1
    index=1;
    for f=fieldnames(data)'
      if index==1, str = [ 'struct(''' f{1} ''',' class2str_eval(data.(f{1}), options) ];
      else         str = [ str  ', ''' f{1} ''', ' class2str_eval(data.(f{1}), options) ]; end
      index=index+1;
    end
    str = [ str ')' ];
  elseif numel(data) > 1
    for index=1:numel(data)
      if index == 1, str = [     '[ ' class2str_eval(data(index), options) ];
      else           str = [ str ', ' class2str_eval(data(index), options) ]; end
    end
    str = [ str ' ]' ];
    if size(data,1) > 1 && size(data,2) == 1, str = [ str '''' ]; end  % transpose
  end
% end class2str_eval

% ------------------------------------------------------------------------------

function str=class2str_validstr(str)
  % validate a string as a single line
  str=strrep(str(:)', sprintf('\n'), ';');
  index = find(str < 32 | str > 127);
  str(index) = ' ';
  str=strrep(str, '''', '''''');