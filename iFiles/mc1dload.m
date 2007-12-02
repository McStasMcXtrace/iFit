function a=mc1dload(a)
% function a=mc1dload(a)
%
% Returns an iData style dataset from a McStas 1d monitor file
%

% Find proper labels for Signal and Axis

xlabel = a.Data.Headers.MetaData.xlabel;
ylabel = a.Data.Headers.MetaData.ylabel;
xlabel(1:length('# xlabel: '))='';
ylabel(1:length('# ylabel: '))='';


Datablock = ['this.' getalias(a,'Signal')];

% First column is the scan parm, we denote that 'x'
setalias(a,'x',[Datablock '(:,1)'],xlabel);
setalias(a,'Signal',[Datablock '(:,2)'],ylabel);
setalias(a,'Error',[Datablock '(:,3)']);
setalias(a,'N',[Datablock '(:,3)'],'# Events');
setaxis(a,1,'x');
