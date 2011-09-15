function a=load_mcstas_1d(a)
% function a=load_mcstas_1d(a)
%
% Returns an iData style dataset from a McStas 1d monitor file, or even simple XYE files
% Some labels are also searched.
%
% Version: $Revision: 1.12 $
% See also: iData/load, iLoad, save, iData/saveas

% inline: load_mcstas_param

% handle input iData arrays
if length(a(:)) > 1
  for index=1:length(a(:))
    a(index) = feval(mfilename, a(index));
  end
  return
end

% Find proper labels for Signal and Axis
a=iData(a);
if isempty(findstr(a,'McStas'))
  warning([ mfilename ': The loaded data set ' a.Tag ' from ' a.Source ' is not a McStas data format.' ]);
  return
end

xlab=''; ylab='';
d = a.Data;
if ~isfield(d,'MetaData'), return; end

if isfield(d,'Headers') && isfield(d.Headers,'MetaData') 
  if ~isempty(findfield(a, 'xlabel')) 
    xlab = a.Data.Headers.MetaData.xlabel;
    xlab(1:max(strfind(xlab,'xlabel')+6))='';
  elseif ~isempty(findfield(a, 'x_label')) 
    xlab = a.Data.Headers.MetaData.x_label;
    xlab(1:max(strfind(xlab,'x_label'))+6)='';
  end

  if ~isempty(findfield(a, 'ylabel')) 
    ylab = a.Data.Headers.MetaData.ylabel;
    ylab(1:max(strfind(ylab,'ylabel')+6))='';
  elseif ~isempty(findfield(a, 'y_label')) 
    ylab = a.Data.Headers.MetaData.y_label;
    ylab(1:max(strfind(ylab,'y_label')+6))='';
  end

  if ~isempty(findfield(a, 'component')) 
    label = strtrim(a.Data.Headers.MetaData.component);
    label(1:length('# component: '))='';
    a.Label = label;
    set(a,'Data.Component', label);
    setalias(a, 'Component', 'Data.Component','Component name');
  end
end

if ~isempty(strfind(a.Title,'McStas 1D monitor'))
  a=xlabel(a, xlab);
  a=title(a, ylab);
end

param = load_mcstas_param(a, 'Param');
a.Data.Parameters = param;
setalias(a, 'Parameters', 'Data.Parameters', 'Instrument parameters');

% ------------------------------------------------------------------------------
% build-up a parameter structure which holds all parameters from the simulation
function param=load_mcstas_param(a, keyword)
  if nargin == 1, keyword='Param:'; end
  param = [];

  par_list = findstr(a, keyword);
  % search strings of the form 'keyword' optional ':', name '=' value
  for index=1:length(par_list)
    line         = par_list{index};
    reversed_line= line(end:-1:1);
    equal_sign_id= find(reversed_line == '=');
    name         = fliplr(strtok(reversed_line((equal_sign_id+1):end),sprintf(' \n\t\r\f;#')));
    if isempty(name)
      column_sign_id = findstr(line, keyword);
      name = strtok(line((column_sign_id+length(keyword)+1):end));
    end
    if isfield(a.Data, name)
      value = getfield(a.Data, name);
    else
      value = strtok(fliplr(reversed_line(1:(equal_sign_id-1))),sprintf(' \n\t\r\f;#'));
      if ~isempty(str2num(value)), value = str2num(value); end
    end
    
    if ~isempty(value) && ~isempty(name) && ischar(name)
      param = setfield(param, name, value);
    end
  end

