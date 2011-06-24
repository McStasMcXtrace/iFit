function s_out = setaxis(a_in,indexes,names,values)
% s = setaxis(s, AxisIndex, AxisName, AxisValues) : set iData axes
%
%   @iData/setaxis function to set iData axes.
%   The function works also when AxisName and AxisIndex are given as cells.
%   The AxisName name must exist in the object (e.g. as alias).
%   When the AxisIndex is empty, the axis is removed, so that
%     setaxis(iData, [], getaxis(iData)) deletes all axis definitions
%   The input iData object is updated if no output argument is specified.
%   If the axis value is specified (as a real value or an alias name), a
%     corresponding Alias is created, and its value is set, so that
%     setaxis(a,1,'x',1:10) creates the 1-st rank axis as the alias 'x'
%     which value is set to 1:10.
%   The Signal/Monitor corresponds to axis 0. Setting its value multiplies 
%     it by the Monitor and then assigns the Signal.
%   Axis 1 is often labeled as 'x' (along rows), 2 as 'y' (along columns), etc...
%   The special syntax a{0} multiplies the value by the Monitor and then assigns 
%   the Signal, and a{n} assigns the axis of rank n.
%     When the assigned value is a char, the axis definition is set.
%     When the assigned value is numeric, the axis value is set (as in set).
%
% input:  s: object or array (iData)
%         AxisIndex: rank of the axis, 
%                    or [] to remove the axis 
%                    or 0 to add a new axis (integer/cell)
%         AxisName: Name of an existing alias/field (char/cellstr)
%         AxisValues: values of the axis (char/alias/numeric)
% output: s: array (iData)
% ex:     setaxis(iData, 1, 'Temperature') defines Temperature as the 'y' axis (rank 1)
%         a{1} =  'Temperature'            does the same
%
% Version: $Revision: 1.17 $
% See also iData, iData/getaxis, iData/get, iData/set, iData/rmaxis

% EF 27/07/00 creation
% EF 23/09/07 iData implementation

s_out=a_in;
if nargin == 1
  % makes a check of axes and Signal, notice invalid ones.
  for index = 1:length(a_in(:))
    a = a_in(index); % current object in array/single element
    for j1=1:length(a.Alias.Axis)
      try
        link = a.Alias.Axis{j1};
        val  = get(a, link);
        if length(find(size(a) > 1)) == 1
          if length(val(:)) ~= size(a, find(size(a) > 1)) & length(val(:)) > 1
            iData_private_warning(mfilename,[ 'the Axis ' a.Alias.Axis{j1} ' ' num2str(j1) '-th rank length ' num2str(length(val(:))) ' does not match the Signal dimension [' num2str(size(a)) '] in object ' inputname(1) ' ' a.Tag '.' ]);
          end
        elseif size(val,j1) ~= size(a, j1) & length(val(:)) ~= size(a, j1)
          iData_private_warning(mfilename,[ 'the Axis ' a.Alias.Axis{j1} ' ' num2str(j1) '-th rank size ' num2str(size(val)) ' does not match the Signal dimension [' num2str(size(a)) '] in object ' inputname(1) ' ' a.Tag '.' ]);
        end
      catch
        iData_private_warning(mfilename,[ 'the Axis ' a.Alias.Axis{j1} ' ' num2str(j1) '-th rank is not valid in object ' inputname(1) ' '  a.Tag '. Defining it.' ]);
        sb.type='{}';
        sb.subs={j1};
        s_out(index)=subsasgn(a, sb, getaxis(a, j1));
      end
    end
  end
  if nargout == 0 & length(inputname(1))
    assignin('caller',inputname(1),s_out);
  end
  return
end
if nargin <= 2
  names=[]; % removes axes
end
if nargin <= 3
  values=[];
end

%if isnumeric(names) & nargin >2, names=mat2str(names(:)); end
if ~isempty(names), 
  if ~isnumeric(names)
    names = cellstr(names); 
  else
    values = names;
  end
else 
  names = { [] }; 
end
if ~iscell(indexes), 
  if isnumeric(indexes), indexes=num2cell(indexes);
  else indexes = { indexes }; end
end
if ~iscell(values), values = { values }; end

s_out = a_in(:);
for i1 = 1:length(s_out)
  a = s_out(i1); % current object in array/single element
  if isnumeric(names) & length(indexes) == 1
    S=struct('type','{}','subs',{indexes}); % change axis as object{rank}=names
    cmd=a.Command;
    a = subsasgn(a, S, names);
    a.Command=cmd;
    a = iData_private_history(a, mfilename, a, indexes, names);
    s_out(i1) = iData(a); % final check
    continue;
  end

  for j1=1:length(names) % loop on axis names
    cmd=a.Command;
    name = names{j1};
    if length(indexes), index = indexes{j1}; else index=[]; end
    
    if ~isempty(name)
      % check that name is not a class member
      f = fieldnames(a);
      if strmatch(lower(name), lower(f))
        iData_private_error(mfilename,[ 'the Alias ' name ' is a protected name in object ' inputname(1) ' ' a.Tag '.' ]);
      end
      
      % does this name already exist as an axis ?
      axis_names = a.Alias.Axis; % this is a cellstr of Axis names
      j2=cellfun('isempty',axis_names);
      if any(j2), axis_names(find(j2)) = ''; end
      
      axis_num   = strmatch(lower(name), lower(axis_names), 'exact');
      
      if length(index)
      if index < 0, index=length(axis_names)+1; end
      end
    else axis_num =[]; end

    % can it be evaluated ?
    try
      val  = get(a, name);
      isvalid=1;
    catch
      isvalid=0;
    end
    
    if isempty(index) & ~isempty(axis_num) % remove existing axis, but do not shift other axis
      a.Alias.Axis(axis_num) = '';
      continue;
    end
    if isempty(index), continue; end
    % assign value to axis when provided
    if nargin==4 & length(name) & ~isempty(index)
      % special case for Signal, which should take into account the Monitor
      if strcmp(name, 'Signal')
        m  = get(a,'Monitor'); m=real(m);
        if not(all(m == 1 | m == 0))
          values{j1} = genop(@times,values{j1},m);
        end
        set(a, 'Signal', values{j1});
      else
        a = setalias(a, name, values{j1});
        %a = setaxis( a, index, name);
      end
      isvalid=1;
    end
    
    if index <= length(a.Alias.Axis) & isempty(name)
      a.Alias.Axis(index) = '';
      continue;
    end

    % axis must exist as an alias/reference to iData field/expression
    if ~isvalid & ~isempty(index)
      todisp='the ';
      if isempty(axis_num)
        todisp = [ todisp 'new ' num2str(index) '-th rank' ];
      else
        todisp = [ todisp 'existing ' num2str(axis_num) '-th rank' ];
      end
      
      a = setalias(a, [ 'Axis_' num2str(index) ], name);
      a = setaxis(a, index,  [ 'Axis_' num2str(index) ]);
      if length(name) > 20, name = [ name(1:17) '...' ];  end
      iData_private_warning(mfilename,[ todisp ' Axis ' name ' is not defined/valid in object ' inputname(1) ' '  a.Tag '.\n\tDefining it as Axis_' num2str(index) ]);
      continue;
    else
    
      if index == 0
        if ~strcmp(name, 'Signal')
          % redefine Signal
          a = setalias(a, 'Signal', name);
        end
      else
        if index <= length(a.Alias.Axis)
          if ~isempty(a.Alias.Axis{index})
            iData_private_warning(mfilename,[ 'redefining Axis ' a.Alias.Axis{index} ' ' num2str(index) '-th rank in object ' inputname(1) ' '  a.Tag ]);
          end
        end
        a.Alias.Axis{index} = name;
      end
    end
    a.Command=cmd;
    a = iData_private_history(a, mfilename, a, index, name);

  end % for alias names

  s_out(i1) = iData(a); % final check
end % for index

if length(s_out) > 1
  s_out = reshape(s_out,size(a_in));
end

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),s_out);
end
