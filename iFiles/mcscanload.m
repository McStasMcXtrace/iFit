function a=mcscanload(a0)
% function a=mcscanload(a0)
%
% Returns iData style datasets from a McStas scan output file
%

% Define alias for the 'raw' datablock
setalias(a0,'Datablock',['this.' getalias(a0,'Signal')]);

% First column is always the scan parm, we denote that 'x'
setalias(a0,'x',['this.' getalias(a0,'Signal') '(:,1)'])
setaxis(a0,1,'x');

siz = size(a0.Signal);
siz = (siz(2)-1)/2;

a = [];
for j=1:siz
  b = a0;
  setalias(b,'Signal',['this.' getalias(a0,'Signal') '(:,' num2str(2*j) ')']);
  setalias(b,'Error',['this.' getalias(a0,'Signal') '(:,' num2str(1+2*j) ')']);
  a = [a b];
end


