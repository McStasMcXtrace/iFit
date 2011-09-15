function a=load_mcstas_2d(a)
% function a=load_mcstas_2d(a)
%
% Returns an iData style dataset from a McStas 2d monitor file
%
% Version: $Revision: 1.8 $
% See also: iData/load, iLoad, save, iData/saveas

% inline: load_mcstas_param

% handle input iData arrays
if length(a(:)) > 1
  for index=1:length(a(:))
    a(index) = feval(mfilename, a(index));
  end
  return
end

a=iData(a);
if isempty(findstr(a,'McStas'))
  warning([ mfilename ': The loaded data set ' a.Tag ' from ' a.Source ' is not a McStas data format.' ]);
  return
end

xlabel=''; ylabel=''; zlabel='';
d = a.Data;
if ~isfield(d,'MetaData'), return; end

if isfield(d,'Headers') && isfield(d.Headers,'MetaData') 
  % Find proper labels for Signal and Axis
  if ~isempty(findfield(a, 'xlabel'))
    xlabel = a.Data.Headers.MetaData.xlabel; 
    xlabel(1:length('# xlabel: '))='';
  else xlabel=''; end
  if ~isempty(findfield(a, 'ylabel'))
    ylabel = a.Data.Headers.MetaData.ylabel;
    ylabel(1:length('# ylabel: '))='';
  else ylabel=''; end
  if ~isempty(findfield(a, 'zlabel'))
    zlabel = a.Data.Headers.MetaData.zlabel;
    zlabel(1:length('# zlabel: '))='';
  else zlabel=''; end
  if ~isempty(findfield(a, 'component')) 
    label = strtrim(a.Data.Headers.MetaData.component);
    label(1:length('# component: '))='';
    a.Label = label;
    set(a,'Data.Component', label);
    setalias(a, 'Component', 'Data.Component','Component name');
  end
end

% Get sizes of x- and y- axes:
siz = size(a.Data.MetaData.variables');
lims = a.Data.MetaData.xylimits;
xax = linspace(lims(1),lims(2),siz(1));
yax = linspace(lims(3),lims(4),siz(2));

% transpose I,E,N fields
if isnumeric(d.MetaData.variables),
  d.MetaData.variables = transpose(d.MetaData.variables);
elseif ischar(d.MetaData.variables)
  link = d.MetaData.variables;
  val  = get(a, link);
  set(a, link, val');
  d = a.Data;
end

% First column is the scan parm, we denote that 'x'
setalias(a,'x',xax,xlabel);
setalias(a,'y',yax,ylabel);
setalias(a,'Signal','Data.MetaData.variables',zlabel);
setalias(a,'I','Signal');
if ~isempty(findfield(a, 'Error'))
  d.MetaData.Errors = transpose(d.MetaData.Errors);
  setalias(a,'Error','Data.MetaData.Errors');
else setalias(a,'Error',0);
end
setalias(a,'E','Error');
if ~isempty(findfield(a, 'Error')) 
  d.MetaData.Events = transpose(d.MetaData.Events);
  setalias(a,'N','Data.MetaData.Events');
end
setaxis(a,1,'x');
setaxis(a,2,'y');

param        = load_mcstas_param(a, 'Param');
d.Parameters = param;
a.Data       = d;
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
