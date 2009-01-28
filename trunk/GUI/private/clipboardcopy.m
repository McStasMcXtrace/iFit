function clipboardData = clipboardcopy(data)
    %clipboardcopy v0.1.0
    %   Usage:
    %       clipboardData = clipboardcopy(data)
    %
    %   Input Arguments:
    %       data (optional)
    %           Data to copy to system clipboard. Data may be 
    %             - a string or an array/cell of strings
    %             - a structure, which is then formated into a Matlab script
    %             - an array of numerics (vector, matrix)
    %             - an image (cdata, numerical array [MxNx3]
    %             - any Matlab object, converted into a structure
    %   Output Arguments:
    %       clipboardData
    %           Struct representing the data on the clipboard. See clipboardpaste for details.
    %   Description:
    %       Sends data to the clipboard. 
    %
    %   Example:
    %       Copy something on to the clipboard. Say that you copy
    %       the number 14584 onto the clipboard.
    %       >> clipboardcopy(14584)
    %       >> clipData = clipboardpaste %now paste it into the workspace
    %       clipData =
    %
    %                  data: '14584'
    %           primaryType: 'text'
    %               subType: 'plain-text'
    %       >> data = clipData.data; %get the actual data from the clipboard
    %
    %   =======================
    %   Written by Bryan Raines on May 5, 2008
    %   Written by E. Farhi on Jan 12 2009 to handle java.awt.images as generated from im2java
    %
    %   See also clipboard, clipboardpaste
    
    usesXClipboard = isunix && ~strcmpi(computer, 'mac');

    if usesXClipboard
        if isempty(get(0,'DefaultFigureXDisplay'))
            error('MATLAB:clipboard:NoXDisplay', 'There is no X display set.');
        end
    else
        err = javachk('awt', 'Clipboard access');
        if (~isempty(err))
            error('MATLAB:clipboard:UnsupportedPlatform', err.message);
        end
    end

    clipboardData     = [];
    error(nargchk(1,1,nargin));
    
    if isempty(inputname(1)), in = [ class(data) '_copy' ];
    else in = inputname(1); end
    
    if ischar(data)
      % COPY: saves as strings
      % 'application/x-java-serialized-object; class=java.lang.String'
      clipboard('copy', data);
      
      % adds \n at all end of lines and make it a single line of chars
      data = cellstr(data);
      data = strcat(data, '<EOL>');
      data = char(data);
      data = transpose(data);
      data = data(:);
      data = transpose(data);
      data = strrep(data, '<EOL>', sprintf('\n'));
      
      ss = java.awt.datatransfer.StringSelection(data); % make the string a Java String
      java.awt.Toolkit.getDefaultToolkit.getSystemClipboard.setContents(ss, []); % copy to clipboard

      clipboardData = clipboardpaste;
    elseif iscellstr(data)
      data = char(data);
      clipboardData = clipboardcopy(data);
    elseif ishandle(data)
      % convert Handle Graphics into an image
      f = getframe(gcf);
      clipboardData = clipboardcopy(f.cdata);
    elseif isnumeric(data)
      sz = size(data);
      if length(sz) == 3
        if sz(3) == 3
          % COPY: makes an image copy
          % 'image/x-java-image; class=java.awt.Image'
          % clipboard('copy', class2str(in, data));
          
          % code from editmenufcn.m
          im = im2java(data); % make a java.awt.image
          jm = javax.swing.ImageIcon(im);
          im_obj = jm.getImage;
          is = com.mathworks.hg.util.ImageSelection(im);
          java.awt.Toolkit.getDefaultToolkit.getSystemClipboard.setContents(is,[]);
          
          clipboardData = clipboardpaste;
          return
        end
      end
      % convert array as a string
      data = mat2str(data);
      clipboardData = clipboardcopy(data);
    else % for structures and objects
      % convert to M code
      data = class2str(in, data);
      clipboardData = clipboardcopy(data);
    end
    
    % --------------------------------------------------------------------------
    
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
  str = [ str 'if ~exist(''iData''), return; end' NL ];
  str = [ str this '_s=' this '; ' this ' = rmfield(' this ',' this '.Alias); ' this ' = iData(' this ');' NL ...
         'setalias(' this ', ' this '_s.Alias.Names, ' this '_s.Alias.Values, ' this '_s.Alias.Labels);' NL ... 
         'setaxis('  this ', mat2str(1:length(' this '_s.Alias.Axis)), ' this '_s.Alias.Axis);' NL ...
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

