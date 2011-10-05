function s_out = set(a_in,varargin)
% [s,...] = set(s, 'PropertyName', Propertyvalue, ...) : set iData properties
%
%   @iData/set function to set iData properties.
%   set(s, 'PropertyName', Propertyvalue, ...}) 
%     sets values into given property names for the iData object s.
%   set(s, Struct.Field, ...)
%     sets values from given structure fields into the iData object s.
%   set(s, CellNames, CellValues, ...)
%     sets values from given cells into the iData object s.
%   set(s) indicates the signification of the iData base properties
%   The input iData object is updated if no output argument is specified.
%
% ex      : set(iData,'Title','A nice Title')
%
% Version: $Revision: 1.9 $
% See also iData, iData/get, iData/setalias, iData/setaxis

% EF 27/07/00 creation
% EF 23/09/07 iData implementation

% calls: setalias

if nargin == 1
  disp('iData object properties:');
  disp('Title:   (string)   title of the Data set');
  disp('Tag:     (string)   unique ID for the Data set');
  disp('Source:  (string)   origin of data (filename/path)');
  disp('Command: (cellstr)  history of commands applied to object');
  disp('Date:    (string)   Data set creation date');
  disp('UserData:(any type) user data storage area');
  disp('Label:   (string)   user label');
  disp('Creator: (string)   application that created this data set');
  disp('User:    (string)   user of this Data set');
  disp('Data:    (any type) Data storage area');
  disp('Axis:    list of axis defined for data set math operations/plotting. Use setaxis/getaxis');
  disp('Alias:   list of aliases/links to data items.                        Use setalias/getalias');
  disp('Signal:  (double)   The signal to be used for data set math operations/plotting.');
  disp('Error:   (double)   The error bars on the signal to be used for data set math operations/plotting.');
  disp('Monitor: (double)   The monitor(statistical weight) on the signal to be used for data set math operations');
  disp('ModificationDate: (string)   last object modification date');
  s_out = iData(a_in);
  return
end

s_out = a_in(:);
for index = 1:length(s_out)
  a = s_out(index); % current object in array/single element
  cmd=a.Command;
  
  i1 = 1; % index in input parameters varargin
  
  while i1<=length(varargin)     % first parse fields and values
    if ischar(varargin{i1})      % normal 'PropertyName', Propertyvalue
      prop_names  = varargin(i1);         % get single PropertyName
      prop_values = varargin(i1+1);       % get single Propertyvalue
      i1 = i1+2;
    elseif isstruct(varargin{i1})       % import structure
      prop_names = fieldnames(varargin{i1});         % get PropertyNames
      prop_values = struct2cell(varargin{i1});
      i1 = i1+1;
    elseif iscell(varargin{i1}) && i1 < length(varargin) % import from 2 cells
      prop_names = varargin{i1};         % get PropertyNames
      prop_values = varargin{i1+1};      % get Propertyvalues
      i1 = i1+2;
    else
      iData_private_error(mfilename, [ 'cannot handle input Property of type ' class(varargin{i1}) ])
    end
    for j1=1:length(prop_names) % loop on properties cell
      prop_name  = prop_names{j1};
      prop_value = prop_values{j1};
      % check for aliases
      alias_names = a.Alias.Names; % this is a cellstr of Alias names
      alias_num   = find(strcmpi(prop_name, alias_names));
      if length(alias_num) >= 1
        alias_values = a.Alias.Values;
        name = alias_names{alias_num(1)};
        prop_name = alias_values{alias_num(1)};
        if ~ischar(prop_name), prop_name=[]; end
        % check if prop_name exists, and use it'
        try
          if (isempty(prop_name) && (isnumeric(prop_value) || islogical(prop_value))) ...
            || isnumeric(prop_name) || strcmp(prop_name, prop_value)
            setalias(a, name, prop_value);
          else
            setalias(a, prop_name, prop_value);
            %eval([ 'a.' prop_name '= prop_value;' ]);
          end
        catch
          if isnumeric(prop_name), 
            prop_name = mat2str(prop_name(1:10)); 
            if length(prop_name) > 15, prop_name = [prop_name(1:12) '...' ]; end 
          end
          iData_private_warning(mfilename, sprintf('can not set Property Alias %s=%s in object %s.', name, prop_name, [ inputname(1) ' '  a.Tag ]))
        end
      else
        if strncmpi(prop_name, 'filename', 8)
          prop_name = 'Source';
        end
        if strncmpi(prop_name, 'alias', 5)
          iData_private_warning(mfilename, 'to set Aliases, use the setalias function');
        elseif strncmpi(prop_name, 'axis', 4)
          iData_private_warning(mfilename, 'to set Axis, use the setaxis function, or the syntax iData{index}.');
        else
          try
            eval(['a.' prop_name '= prop_value;' ]);
          catch
            iData_private_warning(mfilename, sprintf('can not find Property %s in object %s. Setting it.', prop_name, [ inputname(1) ' '  a.Tag ]))
            if ~isempty(prop_value) || isnumeric(prop_value), setalias(a, prop_name, prop_value); end
          end
        end
      end
    end % for properties cell
  end % while
  a.Command=cmd;
  a = iData_private_history(a, mfilename, a, varargin{:});  
  
  s_out(index) = iData(a); % final check
end % for index

if length(s_out) > 1
  s_out = reshape(s_out,size(a_in));
end

if nargout == 0 && ~isempty(inputname(1))
  assignin('caller',inputname(1),s_out);
end
