function a=load_mcstas_1d(a)
% function a=load_mcstas_1d(a)
%
% Returns an iData style dataset from a McStas 1d monitor file
%

% Find proper labels for Signal and Axis

if ~isempty(findfield(a, 'xlabel')) 
  xlabel = a.Data.Headers.MetaData.xlabel;
  xlabel(1:length('# xlabel: '))='';
else xlabel=''; end
if ~isempty(findfield(a, 'ylabel')) 
  ylabel = a.Data.Headers.MetaData.ylabel;
  ylabel(1:length('# ylabel: '))='';
else ylabel=''; end

Datablock = ['this.' getalias(a,'Signal')];

% First column is the scan parm, we denote that 'x'
setalias(a,'x',[Datablock '(:,1)'],xlabel);
setalias(a,'Signal',[Datablock '(:,2)'],ylabel);
if ~isempty(findfield(a, 'Error')) 
  setalias(a,'Error',[Datablock '(:,3)']);
else
  setalias(a,'Error',0);
end
setalias(a,'E','Error');
if ~isempty(findfield(a, 'Error')) 
  setalias(a,'N',[Datablock '(:,3)'],'# Events');
end
setaxis(a,1,'x');

