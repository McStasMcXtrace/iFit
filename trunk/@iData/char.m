function c = char(s)
% c = char(s) : convert iData into character
%
%   @iData/char: function to convert iData objects into char
%   returns the iData title/filename
%
% input:  s: object or array (iData) 
% output: c: iData identification (char)
%  
% See also  iData/cell, iData/double, iData/struct, 
%           iData/char, iData/size
%

% EF 23/09/07 iData implementation

c=[];
for index=1:length(s)
  t = s(index);
  T = t.Title;  if iscell(T), T = T{1}; end
  cmd = t.Command{end};
  if length(cmd) > 23, cmd = [ cmd(1:20) '...' ]; end
  c = strvcat(c, [ 'iData ' t.Tag '=' cmd ' [' num2str(size(t)) '] "' deblank(T) '" <' t.Source '>' ]); 
end

