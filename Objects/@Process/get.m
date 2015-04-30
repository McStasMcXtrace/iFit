function [varargout] = get(this,varargin)
% [...] = get(s, 'PropertyName', ...) : get Process object properties
%
%   @Process/get function to get Process properties.
%   get(s) displays all property names and their current values for
%     the Process object 's'.
%   get(s,'PropertyName',...) returns only particular properties.
%     The PropertyName can be any Process object field, or a model parameter name
%       or 'p' to designate the vector of parameter values (when previously set).
%   Input 's' can be a single Process or a Process array
%
% input:  s: object or array (Process)
%         PropertyName: name of Property to search (char)
% output: property: property value in 's' (cell)
% ex :    get(Process) or get(Process,'command')
%
% Version: $Revision$
% See also Process, Process/disp, Process/refresh, Process/exit

persistent fields

if isempty(fields), fields=fieldnames(Process); end

% handle array of objects
varargout = {};

if numel(this) > 1
  varg = {};
  for index=1:numel(this)
    varg{end+1} = get(this(index), varargin{:});
  end
  varargout{1} = varg;
  return
end

if nargin == 1
  disp(this, inputname(1));
  varargout{1} = display(this, inputname(1));
  return
end

% handle single object
for index=1:length(varargin)
  property = varargin{index}; % get PropertyName
  if isempty(property), continue; end
  if ~ischar(property)
    error([ mfilename ': PropertyName should be char strings in object ' inputname(1) ' ' this.Tag ' "' this.Title '" and not ' class(property) ]);
  end
  % test if this is a unique property, or a composed one
  if isvarname(property)  % extract Process field/alias
    if any(strcmp(property, fields))
      b = this.(property);               % direct static field
      if isnumeric(b) && any(strcmp(property, {'Date','ModificationDate'}))
        b = datestr(b);
      end
      varargout{1} = b;
    end
  end
end

if isempty(varargout)
  varargout={[]};
end
