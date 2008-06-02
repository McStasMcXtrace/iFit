function a=load_mcstas_scan(a0)
% function a=load_mcstas_scan(a0)
%
% Returns iData style datasets from a McStas scan output file
%
a=iData(a);
% Define alias for the 'raw' datablock
setalias(a0,'Datablock',['this.' getalias(a0,'Signal')]);

if ~isempty(findfield(a0, 'xlabel')) 
  xlabel = a0.Data.Headers.MetaData.xlabel;
  xlabel(1:length('# xlabel: '))='';
else xlabel=''; end
if ~isempty(findfield(a0, 'ylabel')) 
  ylabel = a0.Data.Headers.MetaData.ylabel;
  ylabel(1:length('# ylabel: '))='';
else ylabel=''; end

% get the column labels
cnames=strread(a0.Data.Headers.MetaData.variables,'%s','delimiter',' ');
cnames=cnames(3:end);
xlabel=cnames(1);

% First column is always the scan parm, we denote that 'x'
setalias(a0,'x',['this.' getalias(a0,'Signal') '(:,1)'],xlabel);
setaxis(a0,1,'x');

siz = size(a0.Signal);
siz = (siz(2)-1)/2;

a = [];
for j=1:siz
  b = copyobj(a0);
  ylabel=cnames(2*j);
  setalias(b,'Signal',['this.' getalias(a0,'Signal') '(:,' num2str(2*j) ')'],ylabel);
  if ~isempty(findfield(a0, '_ERR')) 
    setalias(b,'Error',['this.' getalias(a0,'Signal') '(:,' num2str(1+2*j) ')']);
  end
  b.Title=[ char(ylabel) ': ' char(b.Title) ];
  a = [a b];
end


