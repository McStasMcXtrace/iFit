function a=load_mcstas_scan(a0)
% function a=load_mcstas_scan(a0)
%
% Returns iData style datasets from a McStas scan output file
%
a=iData(a0);
% Define alias for the 'raw' datablock
setalias(a0,'Datablock',['this.' getalias(a0,'Signal')]);

% get the column labels
cnames=strread(a0.Data.Headers.MetaData.variables,'%s','delimiter',' ');
cnames=cnames(3:end);

if ~isempty(findfield(a0, 'xlabel')) 
  xlabel = deblank(a0.Data.Headers.MetaData.xlabel);
  xlabel(1:length('# xlabel: '))='';
else xlabel=''; end
if ~isempty(findfield(a0, 'ylabel')) 
  ylabel = deblank(a0.Data.Headers.MetaData.ylabel);
  ylabel(1:length('# ylabel: '))='';
else ylabel=''; end
if ~isempty(findfield(a0, 'xvars')) 
  xvars = deblank(a0.Data.Headers.MetaData.xvars);
  xvars(1:length('# xvars: '))='';
else xvars=''; end

if ~isempty(xvars)
  xvars_i = find(cellfun('isempty', strfind(cnames,xvars)) == 0);
  if ~isempty(xvars_i)
    setalias(a0,'x',['this.' getalias(a0,'Signal') '(:,' num2str(xvars_i) ')' ],xvars); % create an alias for xvars
    setalias(a0,xvars,['this.' getalias(a0,'Signal') '(:,' num2str(xvars_i) ')' ],xvars); % create an alias for xvars
    % locate xvars label and column
    xlabel=xvars;
  end

  % Define scanning variable
  setaxis(a0,1,'x');
end

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
  b.Title = [ char(ylabel) ': ' char(b.Title) ];
  b.Label = [ char(ylabel) '(' xvars ')' ];
  a = [a b];
end


