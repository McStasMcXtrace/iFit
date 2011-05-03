function a=load_mcstas_1d(a)
% function a=load_mcstas_1d(a)
%
% Returns an iData style dataset from a McStas 1d monitor file, or even simple XYE files
% Some labels are also searched.
%

% Find proper labels for Signal and Axis
a=iData(a);

if ~isempty(findfield(a, 'xlabel')) 
  xlab = a.Data.Headers.MetaData.xlabel;
  xlab(1:max(strfind(xlab,'xlabel')+6))='';
elseif ~isempty(findfield(a, 'x_label')) 
  xlab = a.Data.Headers.MetaData.x_label;
  xlab(1:max(strfind(xlab,'x_label'))+6)='';
else xlab=''; end

if ~isempty(findfield(a, 'ylabel')) 
  ylab = a.Data.Headers.MetaData.ylabel;
  ylab(1:max(strfind(ylab,'ylabel')+6))='';
elseif ~isempty(findfield(a, 'y_label')) 
  ylab = a.Data.Headers.MetaData.y_label;
  ylab(1:max(strfind(ylab,'y_label')+6))='';
else ylab=''; end

if ~isempty(findfield(a, 'component')) 
  label = strtrim(a.Data.Headers.MetaData.component);
  label(1:length('# component: '))='';
  a.Label = label;
  set(a,'Data.Component', label);
  setalias(a, 'Component', 'Data.Component','Component name');
end

% special case for McStas files and XYE (2-4 columns) files
n = size(a,2); % number of columns
if (n >= 2 && n <= 4 && size(a,1) >= 5) || ~isempty(strfind(a.Title,'McStas 1D monitor'))

  Datablock = ['this.' getalias(a,'Signal')];

  % First column is the scan parm, we denote that 'x'
  setalias(a,'x',[Datablock '(:,1)'],xlab);
  setalias(a,'Signal',[Datablock '(:,2)'],ylab);
  if n>=3
    setalias(a,'Error',[Datablock '(:,3)']);
  else
    setalias(a,'Error',[]);
  end
  setalias(a,'E','Error');
  if ~isempty(findfield(a, 'Error')) || n >= 4
    setalias(a,'N',[Datablock '(:,4)'],'# Events');
  end
  setaxis(a,1,'x');
else
  xlabel(a, xlab);
  ylabel(a, ylab);
end

