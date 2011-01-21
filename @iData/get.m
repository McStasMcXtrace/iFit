function [varargout] = get(a_in,varargin)
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
% Version: $Revision: 1.15 $
% See also iData, iData/set, iData/getalias, iData/getaxis, iData/findobj

% EF 27/07/00 creation
% EF 23/09/07 iData implementation
% ============================================================================
% private function: iData_private_getalias

for index = 1:length(a_in(:)) % works with object arrays
  argout = 1;
  s = a_in(index);
  s_cell=struct2cell(struct(s));  
  names = fieldnames(struct(s));

  if nargin == 1   
    if length(a_in(:)) == 1
      disp(s);
    end
    varout{argout,index} = s; % the display method is used for get(object)
    argout = argout + 1;
  else
    for i=1:length(varargin)
      val = '';
      if ischar(varargin{i})
        fieldname = varargin{i}; % get PropertyName
      else
        iData_private_error(mfilename, [ 'PropertyNames should be char strings in object ' inputname(1) ' ' s.Tag ' and not ' class(varargin{i}) ]);
      end
      % test private/protected fields
      if strcmp(lower(fieldname), 'alias') | strncmp(lower(fieldname), 'alias.', 6)
        varout{argout,index} = getalias(s); % the display method is used for get(iData)
        argout = argout + 1;
        % warning('iData/get: the Alias field name is protected. Use setalias/getalias to handle Aliases.');
        continue;
      end
      if strcmp(lower(fieldname), 'axis') | strncmp(lower(fieldname), 'axis.', 5)
        varout{argout,index} = getaxis(s); % the display method is used for get(iData)
        argout = argout + 1;
        % warning('iData/get: the Axis field name is protected. Use setaxis/getaxis to handle Axis.');
        continue;
      end

      % get property from real fields (not Alias/Axis)
      field_not_found=1;
      if isempty(fieldname), val=s; field_not_found=0; end
      if field_not_found  % searches within all aliases
        link = ''; name = '';
        % if not found, searches for user Aliases
        if field_not_found
          alias_names = s.Alias.Names; % this is a cellstr of Alias names
          alias_num   = strmatch(lower(fieldname), lower(alias_names), 'exact');
          if ~isempty(alias_num)
            alias_values = s.Alias.Values;
            name = alias_names{alias_num(1)};
            link = alias_values{alias_num(1)};
            field_not_found=0;
            if ~isempty(s.Data)
              val = iData_getalias(s, link, name); % eval link value
            else
              val = [];
            end
          end
        end
      end
      if field_not_found  % searches within iData fields
        cell_num = strmatch(lower(fieldname), lower(names), 'exact');
        if ~isempty(cell_num)        % was the property found immediately (not a sub-field) ?
          val = s_cell{cell_num(1)}; % get field from structure
          field_not_found=0;
        else
          try  % try eval (which works with sub-fields, but is case sensitive)
            val = eval([ 's.' fieldname ';' ]);
            field_not_found=0;
          catch
          end
        end
      end
      if field_not_found
        iData_private_error(mfilename, sprintf('can not find Property "%s" in object %s.', fieldname, [ inputname(1) ' ' s.Tag ] ));
        val=[];
      end
      varout{argout,index} = val;
      argout = argout + 1;
    end % for varargin
  end % if nargin = 1
end % for index

if ~length(a_in(:)) == 1
  varargout{1}=[];
  return
end
for i = 1:(argout-1)
  if length(a_in(:)) == 1
    varargout{i} = varout{i};
  else
    varargout{i} = reshape({varout{i,:}},size(a_in));
  end
end

% ============================================================================
% private function iData_getalias
function val = iData_getalias(this,link,name)
% iData_getalias: iData alias evaluation
%   evaluates s.name to be first s.link, then 'link' (with 'this' defined).
%   NOTE: for standard Aliases (Error, Monitor), makes a dimension check on Signal

% EF 23/09/07 iData impementation

val = [];
if (isnumeric(link) | islogical(link)) & ~isempty(link), 
  if strcmp(name, 'Monitor') && all(val == 0)
    val = 1;
  else
    val = link; 
  end
  return; 
end
if strcmp(link, name), return; end       % avoids endless iteration.
if ~isempty(link)
  try
    val = get(this, link);               % link may refer to other Aliases
  catch
    try
      val = eval([ '[ ' link ' ]' ]);
    catch
      lasterr
      iData_private_error(mfilename,[ 'can not evaluate Alias "' name '" as ''' link ''' in iData object ' this.Tag ]);
    end
  end
end

% link value has been evaluated, do check in case of standard aliases
if strcmp(name, 'Error')  % Error is sqrt(Signal) if not defined
  if isempty(val) & isnumeric(get(this,'Signal'))
    try
    val = sqrt(abs(get(this,'Signal')));
    catch
    val=0;
    end
  end
  if length(val) ~= 1 & ~all(size(val) == size(this))
    iData_private_warning(mfilename,[ 'The Error [' num2str(size(val)) '] has not the same size as the Signal [' num2str(size(this)) '] in iData object ' this.Tag '.\n\tTo use the default Error=sqrt(Signal) use s.Error=[].' ]);
  end
elseif strcmp(name, 'Monitor')  % monitor is 1 by default
  if isempty(val) | all(val == 0)
    val = ones(size(this));
  end
  if length(val) ~= 1 & ~all(size(val) == size(this))
    iData_private_warning(mfilename,[ 'The Monitor [' num2str(size(val)) '] has not the same size as the Signal [' num2str(size(this)) '] in iData object ' this.Tag '.\n\tTo use the default Monitor=1 use s.Monitor=[].' ]);
  end
end

