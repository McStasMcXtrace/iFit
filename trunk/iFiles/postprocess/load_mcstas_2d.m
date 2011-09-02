function a=load_mcstas_2d(a)
% function a=load_mcstas_2d(a)
%
% Returns an iData style dataset from a McStas 2d monitor file
%
a=iData(a);
if isempty(findstr(a,'McStas'))
  warning([ mfilename ': The loaded data set ' a.Tag ' is not a McStas data format.' ]);
  return
end

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

% Get sizes of x- and y- axes:
siz = size(a.Data.MetaData.variables');
lims = a.Data.MetaData.xylimits;
xax = linspace(lims(1),lims(2),siz(1));
yax = linspace(lims(3),lims(4),siz(2));

% First column is the scan parm, we denote that 'x'
setalias(a,'x',xax,xlabel);
setalias(a,'y',yax,ylabel);
setalias(a,'Signal','transpose(this.Data.MetaData.variables)',zlabel);
setalias(a,'I','Signal');
if ~isempty(findfield(a, 'Error')) 
  setalias(a,'Error','transpose(this.Data.MetaData.Errors)');
else setalias(a,'Error',0);
end
setalias(a,'E','Error');
if ~isempty(findfield(a, 'Error')) 
  setalias(a,'N','transpose(this.Data.MetaData.Events)');
end
setaxis(a,1,'x');
setaxis(a,2,'y');

