function a=load_mcstas_1d(a)
% function a=load_mcstas_1d(a)
%
% Returns an iData style dataset from a McStas 1d monitor file
%

% Find proper labels for Signal and Axis
a=iData(a);

if ~isempty(findfield(a, 'xlabel')) 
  xlabel = a.Data.Headers.MetaData.xlabel;
  xlabel(1:length('# xlabel: '))='';
else xlabel=''; end
if ~isempty(findfield(a, 'ylabel')) 
  ylabel = a.Data.Headers.MetaData.ylabel;
  ylabel(1:length('# ylabel: '))='';
else ylabel=''; end
if ~isempty(findfield(a, 'component')) 
  label = strtrim(a.Data.Headers.MetaData.component);
  label(1:length('# component: '))='';
  a.Label = label;
  set(a,'Data.Component', label);
  setalias(a, 'Component', 'Data.Component','Component name');
end

Datablock = ['this.' getalias(a,'Signal')];

% First column is the scan parm, we denote that 'x'
setalias(a,'x',[Datablock '(:,1)'],xlabel);
setalias(a,'Signal',[Datablock '(:,2)'],ylabel);
try
  setalias(a,'Error',[Datablock '(:,3)']);
catch
  setalias(a,'Error',0);
end
setalias(a,'E','Error');
if ~isempty(findfield(a, 'Error')) 
  setalias(a,'N',[Datablock '(:,3)'],'# Events');
end
setaxis(a,1,'x');

