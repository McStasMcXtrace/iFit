function [varargout] = get(this,varargin)
% [...] = get(s, 'PropertyName', ...) : get iData object properties
%
%   @iData/get function to get iData properties.
%   get(s) displays all property names and their current values for
%     the iData object 's'.
%   get(s,'PropertyName',...) returns only particular properties.
%     the PropertyName may also be an object Alias or an Axis
%   Input 's' can be a single iData or a iData array
%
% input:  s: object or array (iData)
%         PropertyName: name of Property to search (char)
% output: property: property value in 's' (cell)
% ex :    get(iData) or get(iData,'Title')
%
% Version: $Revision: 1.27 $
% See also iData, iData/set, iData/getalias, iData/getaxis, iData/findobj

% EF 27/07/00 creation
% EF 23/09/07 iData implementation
% ============================================================================
% calls: subsref

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
  disp(this);
  varargout{1} = display(this);
  return
end

% handle single object
for index=1:length(varargin)
  property = varargin{index}; % get PropertyName
  if isempty(property), continue; end
  if ~ischar(property)
    iData_private_error(mfilename, [ 'PropertyName should be char strings in object ' inputname(1) ' ' this.Tag ' and not ' class(property) ]);
  end
  % test if this is a unique property, or a composed one
  if isvarname(property)  % extract iData field/alias
    s = struct('type','.','subs', property);
    varargout{1} = subsref(this, s);              % calls subsref directly (single subsref level)
  else % this is a compound property, such as get(this,'Data.Signal')
    try
      varargout{1} = eval([ 'this.' property ]);  % calls subsref by eval (recusive subsref levels)
    catch
      varargout{1} = eval([ property ]);          % this is a full expression: evaluate it...
    end
  end
end

