function t = title(a, name)
% b = title(s) : get/set the model Name
%
%   @iFunc/title function to get and set the model Name/Title
%
%   t = title(s)
%      get the current model name/title
%   t = title(s, t)
%      set the model Name to 't'
%
% input:  s: object or array (iFunc)
% output: t: model named (char)
% ex:     b=title(gauss);
%
% Version: $Date$
% See also iFunc

if nargin == 1
  t = a.Name;
elseif ischar(name) || iscellstr(name)
  t = char(name);
  a.Name = t;
  
  if nargout == 0 && ~isempty(inputname(1)) && isa(a,'iFunc')
    assignin('caller',inputname(1),a);
  end
end

