function a=load_mcstas_2d(a)
% function a=load_mcstas_2d(a)
%
% Returns an iData style dataset from a McStas 2d monitor file
%

% Find proper labels for Signal and Axis

xlabel = a.Data.Headers.MetaData.xlabel;
ylabel = a.Data.Headers.MetaData.ylabel;
zlabel = a.Data.Headers.MetaData.zlabel;
xlabel(1:length('# xlabel: '))='';
ylabel(1:length('# ylabel: '))='';
zlabel(1:length('# zlabel: '))='';

% Get sizes of x- and y- axes:
siz = size(a.Data.MetaData.variables);
lims = a.Data.MetaData.xylimits;
xax = linspace(lims(1),lims(2),siz(1));
yax = linspace(lims(3),lims(4),siz(2));

% First column is the scan parm, we denote that 'x'
setalias(a,'x',xax,xlabel);
setalias(a,'y',yax,ylabel);
setalias(a,'Signal',a.Data.MetaData.variables,zlabel);
setalias(a,'I','Signal');
setalias(a,'Error',a.Data.MetaData.Errors);
setalias(a,'E','Error');
setalias(a,'N',a.Data.MetaData.Events);
setaxis(a,1,'x');
setaxis(a,2,'y');
