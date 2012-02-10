function a=load_mcstas_1d(a)
% function a=load_mcstas_1d(a)
%
% Returns an iData style dataset from a McStas 1d/2d/list monitor file
% as well as simple XYE files
% Some labels are also searched.
%
% Version: $Revision: 1.14 $
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
    xlab(1:max(strfind(xlab,'xlab')+6))='';
  elseif ~isempty(findfield(a, 'x_label')) 
    xlab = a.Data.Headers.MetaData.x_label;
    xlab(1:max(strfind(xlab,'x_label'))+6)='';
  else xlab='';
  end

  if ~isempty(findfield(a, 'ylabel')) 
    ylab = a.Data.Headers.MetaData.ylabel;
    ylab(1:max(strfind(ylab,'ylab')+6))='';
  elseif ~isempty(findfield(a, 'y_label')) 
    ylab = a.Data.Headers.MetaData.y_label;
    ylab(1:max(strfind(ylab,'y_label')+6))='';
  else ylab='';
  end
  
  if ~isempty(findfield(a, 'zlabel')) 
    zlab = a.Data.Headers.MetaData.zlabel;
    zlab(1:max(strfind(zlab,'zlab')+6))='';
  elseif ~isempty(findfield(a, 'z_label')) 
    zlab = a.Data.Headers.MetaData.z_label;
    zlab(1:max(strfind(zlab,'z_label')+6))='';
  else zlab='';
  end
 
  if ~isempty(findfield(a, 'component')) 
    label = strtrim(a.Data.Headers.MetaData.component);
    label(1:length('# component: '))='';
    a.Label = label;
    a.Data.Component = label;
    setalias(a, 'Component', 'Data.Component','Component name');
  end
  
  if ~isempty(findfield(a, 'Creator'))
    creator = a.Data.Headers.MetaData.Creator;
    creator(1:length('# Creator: '))='';
    a.Creator=creator; 
  end
end

% treat specific data formats 1D, 2D, List for McStas ==========================
if ~isempty(strfind(a.Format,'McStas 1D monitor'))
  xlabel(a, xlab);
  title(a, ylab);
elseif ~isempty(strfind(a.Format,'McStas 2D monitor'))
  % Get sizes of x- and y- axes:
  siz = size(a.Data.MetaData.variables');
  lims = a.Data.MetaData.xylimits;
  xax = linspace(lims(1),lims(2),siz(1));
  yax = linspace(lims(3),lims(4),siz(2));

  % First column is the scan parm, we denote that 'x'
  setalias(a,'y',xax,xlab);
  setalias(a,'x',yax,ylab);
  setalias(a,'Signal','Data.MetaData.variables',zlab);
  setalias(a,'I','Signal');
  if ~isempty(findfield(a, 'Error'))
    setalias(a,'Error','Data.MetaData.Errors');
  else setalias(a,'Error',0);
  end
  setalias(a,'E','Error');
  if ~isempty(findfield(a, 'Error')) 
    setalias(a,'N','Data.MetaData.Events');
  end
  setaxis(a,1,'x');
  setaxis(a,2,'y');
elseif ~isempty(strfind(a.Format,'McStas list monitor'))
  % the Signal has been set to the biggest field, which contains indeed the List
  list = getalias(a, 'Signal');
  setalias(a, 'List', list, 'List of events');

  % column signification is given by tokens from the ylab
  columns = strread(ylab,'%s','delimiter',' ');
  index_axes = 0;
  for index=1:length(columns)
    setalias(a, columns{index}, [ list '(:,' num2str(index) ')' ]);
    if strcmp(columns{index}, 'p')
      setalias(a, 'Signal', columns{index}, 'Intensity');
    elseif index_axes < 3
      index_axes = index_axes +1;
      setaxis(a, index_axes, columns{index});
    end
  end
end

% get the instrument parameters
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

