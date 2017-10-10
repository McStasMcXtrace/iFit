function out = opensim(filename)
%OPENSIM Open a McStas SIM file, display it
%        and set the 'ans' variable to an iData object with its content
% (c) E.Farhi, ILL. License: EUPL.

if ~isa(filename,'iData')
  out = iData(iLoad(filename,mfilename)); % no post-processing
else
  out = filename;
end
clear filename;

if numel(out) > 1
  % handle input iData arrays
  for index=1:numel(out)
    out(index) = feval(mfilename, out(index));
  end
end

if numel(out) == 1
  if ~isempty(findstr(out,'simulation description file'))
    % this is a McCode file
    
    % Find filename fields in sim struct:
    filenames = findstr(out,'filename');
    if ~iscellstr(filenames), filenames = { filenames }; end
    dirname   = fileparts(out.Source);

    a=[];
    if length(filenames(:)) > 0
      % This is a McCode 'overview' plot
      for j=1:length(filenames(:))
        filename = filenames{j};
        filename(~isstrprop(filename,'print')) = '';
        filename(1:length('filename: '))='';
        filename=strtrim(filename);
        filename(findstr(' ',filename):length(filename))='';
        if isempty(filename), continue; end
        filename = fullfile(dirname,filename);
        a = [ a iData(filename) ];
      end
    else
      % This is a .sim from a scan
      filename = 'mccode.dat';
      filename = fullfile(dirname,filename);
      if exist(filename)
        a = iData(filename);
      end
    end
    out = a;
    clear a;
  elseif ~isempty(findstr(out,'Isotropic_Sqw'))
    out = opensqw(out);
  elseif ~isempty(findstr(out, 'multiarray_1d'))
    out = load_mcstas_scan(out);
  else
    out = load_mcstas_dat(out); % private, below
  end
  
  if ~isempty(out) && isempty(findstr(out,'McStas')) && isempty(findstr(out,'McCode'))
    warning([ mfilename ': The loaded data set ' out.Tag ' from ' out.Source ' is not a McCode data format.' ]);
    return
  end

end

if ~nargout
  figure; subplot(out);
  
  if ~isdeployed
    assignin('base','ans',out);
    ans = out
  end
end

% ------------------------------------------------------------------------------
function a=load_mcstas_dat(a)
% function a=load_mcstas_dat(a)
%
% Returns an iData style dataset from a McCode 1d/2d/list monitor file
% as well as simple XYE files
% Some labels are also searched.
%
% Version: $Date$
% See also: iData/load, iLoad, save, iData/saveas

% inline: load_mcstas_param

% Find proper labels for Signal and Axis

xlab=''; ylab=''; zlab='';
d = a.Data;

i = findfield(a, {'xlabel','x_label'}, 'char exact');
if ~isempty(i)
  if iscell(i) l=get(a,i{1}); else l=get(a,i); end
  [l, xlab] = strtok(l, ':='); xlab=strtrim(xlab(2:end));
  if isempty(xlab), xlab=l; end
end

i = findfield(a, {'ylabel','y_label'}, 'char exact cache');
if ~isempty(i)
  if iscell(i) l=get(a,i{1}); else l=get(a,i); end
  [l, ylab] = strtok(l, ':='); ylab=strtrim(ylab(2:end));
  if isempty(ylab), ylab=l; end
end

i = findfield(a, {'xlabel','x_label'}, 'char exact cache');
if ~isempty(i)
  if iscell(i) l=get(a,i{1}); else l=get(a,i); end
  [l, zlab] = strtok(l, ':='); zlab=strtrim(zlab(2:end));
  if isempty(zlab), zlab=l; end
end

i = findfield(a, 'component','char exact cache');
if ~isempty(i)
  if iscell(i) l=get(a,i{1}); else l=get(a,i); end
  [l, label] = strtok(l, ':='); label=strtrim(label(2:end));
  if isempty(label), label=l; end
  a.Label = label;
  a.Data.Component = label;
  setalias(a, 'Component', 'Data.Component','Component name');
end

i = findfield(a, 'Creator','char exact cache'); % should have at least the iData.Creator property
if iscell(i) && length(i) > 1
  l=get(a,i{end});
  [l, creator] = strtok(l, ':='); creator=strtrim(creator(2:end));
  if isempty(creator), creator=l; end
  a.Creator=creator; 
end

i = findfield(a, 'type','char exact cache'); % should have at least the iData.Creator property
if ~isempty(i)
  if iscell(i) t=get(a,i{1}); else t=get(a,i); end
  a.Format = [ a.Format ' ' t];
end

clear d

% check that guessed Signal is indeed what we look for
signal = getalias(a, 'Signal');
if ischar(signal) && (~isempty(strfind(signal, 'MetaData')) ...
    || ~isempty(strfind(signal, 'errors')) ...
    || ~isempty(strfind(signal, 'events')))
  % biggest field is not the list but some MetaData, search other List 
  % should be 'Data.MataData.variables' or 'Data.I'
  for search = {'Data.I','Data.Sqw','Data.MetaData.variables','Data.data'}
    if ~isempty(findfield(a,search,'case exact numeric cache')), signal = search; break; end
  end
  if isempty(signal)
    [match, types, dims] = findfield(s, '', 'numeric cache');
    if length(match) > 1, signal = match{2}; end
  end
  if ~isempty(signal), setalias(a, 'Signal', signal); end
end
siz = size(getaxis(a,'Signal')');

% treat specific data formats 1D, 2D, List for McStas ==========================
if ~isempty(strfind(a.Format,'0D monitor'))
  a = setalias(a, 'Signal', 'Data.values(1)');
  a = setalias(a, 'Error' , 'Data.values(2)');
  a = setalias(a, 'I', 'Signal');
  a = setalias(a, 'E', 'Error');
  a = setalias(a, 'N', 'Data.values(3)');
elseif (~isempty(strfind(a.Format,'1D monitor')) || ~isempty(strfind(a,'array_1d'))) && size(a,2) ~= 4
  xlabel(a, xlab);
  title(a, ylab);
  a = setalias(a, 'I', 'Signal');
  a = setalias(a, 'E', 'Error');
  siz = numel(getaxis(a,'Signal')');
  % set axes
  i = findfield(a, {'xylimits','xlimits'}, 'numeric exact cache');
  if iscell(i) && ~isempty(i), i = i{1}; end
  lims = get(a,i);
  if numel(lims)>= 2, 
    xax = linspace(lims(1),lims(2),siz);
    setalias(a,'x',xax,xlab);
    setaxis(a,1,'x');
  end
elseif ~isempty(strfind(a,'array_1d')) && size(a,2) == 4
  % recover a wrongly imported McCode 1D, using header, and set default columns
  sig = getalias(a,'Signal');
  setalias(a,'Axis_1', [ sig '(:,1)' ]);
  setalias(a,'Signal', [ sig '(:,2)' ]);
  setalias(a,'Error',  [ sig '(:,3)' ]);
  setalias(a,'Events', [ sig '(:,4)' ]);
  setaxis(a,1,'Axis_1');
  header = findfield(a,'Attributes','char');
  if ~isempty(header)
    if iscell(header), header = get(a, header{1});
    else               header = get(a, header); end
    header = textscan(header, '%s','Delimiter','#');
    header = str2struct(header{1});
    setalias(a, 'Header', header);
    a.Format = 'McCode 1D monitor';
    a = opensim(a); % this will now make a better job
    return
  end
elseif (~isempty(strfind(a.Format,'2D monitor')) || ~isempty(strfind(a,'array_2d'))) ...
    && isempty(strfind(a.Format,'list monitor'))
  % Get sizes of x- and y- axes:
  i = findfield(a, 'variables', 'numeric exact cache');
  if iscell(i) && ~isempty(i), i = i{1}; end
  if numel(get(a, i)) >= prod(siz)    % e.g. i='Data.MetaData.variables'
    setalias(a,'Signal',i,zlab);
  end
  siz = size(getaxis(a,'Signal')');
  setalias(a,'I','Signal');
  i = findfield(a, 'Errors', 'numeric exact cache');
  if ~isempty(i)
    setalias(a,'Error',i{1});
  else setalias(a,'Error',0);
  end
  setalias(a,'E','Error');
  i = findfield(a, 'Events', 'numeric exact cache');
  if ~isempty(i) 
    setalias(a,'N',i{1});
  end
  
  % set axes
  i = findfield(a, {'xylimits','xlimits'}, 'numeric exact cache');
  if iscell(i) && ~isempty(i), i = i{1}; end
  lims = get(a,i);
  if numel(lims)>= 2, 
    xax = linspace(lims(1),lims(2),siz(1));
    setalias(a,'y',xax,xlab);
    setaxis(a,2,'y');
  end
  if numel(lims)>= 4, 
    yax = linspace(lims(3),lims(4),siz(2));
    setalias(a,'x',yax,ylab);
    setaxis(a,1,'x');
  end

elseif ~isempty(strfind(a.Format,'list monitor'))
  % the Signal should contain the List
  list = getalias(a, 'Signal');
  if ischar(list)
    setalias(a, 'List', list, 'List of events');

    % column signification is given by tokens from the ylab
    columns = strread(ylab,'%s','delimiter',' ');
    index_axes = 0;
    for index=1:length(columns)
      setalias(a, columns{index}, [ list '(:,' num2str(index) ')' ]);
      if index==1
        setalias(a, 'Signal', columns{index});
      elseif index_axes < 3
        index_axes = index_axes +1;
        setaxis(a, index_axes, columns{index});
      end
    end
    if ~isfield(a, 'N'), setalias(a, 'N', length(a{0})); end
  end
end

% build the title: 
%   sum(I) sqrt(sum(I_err^2)) sum(N)


values = get(a, findfield(a, 'values','cache'));
if ~isempty(values)
  if iscell(values) values=values{1}; end
  t_sum = sprintf(' I=%g I_err=%g N=%g', values);
else
  t_sum = ''; 
end
%   X0 dX, Y0 dY ...

t_XdX = '';
if ~isscalar(a)
  if ndims(a) == 1; ax='X'; else ax = 'YXZ'; end
  for index=1:ndims(a)
    [dx,x0]=std(a,index);
    t_XdX = [t_XdX sprintf(' %c0=%g d%c=%g;', ax(index), x0, ax(index), dx) ];
  end
end
a.Title = [ a.Title, t_sum, t_XdX ];
setalias(a,'statistics',  t_XdX,'Center and Gaussian half width');
setalias(a,'values',values,'I I_err N');

% get the instrument parameters
param = load_mcstas_param(a, 'Param');
a.Data.Parameters = param;
setalias(a, 'Parameters', 'Data.Parameters', 'Instrument parameters');

% end of loader

% ------------------------------------------------------------------------------
function a=load_mcstas_scan(a0)
% function a=load_mcstas_scan(a0)
%
% Returns iData style datasets from a McStas scan output file
%
% See also: iData/load, iLoad, save, iData/saveas

% Define alias for the 'raw' datablock
setalias(a0,'Datablock',['this.' getalias(a0,'Signal')]);

% get the column labels
cnames=strread(a0.Data.Attributes.MetaData.variables,'%s','delimiter',' ');
cnames=cnames(3:end);

if ~isempty(findfield(a0, 'xlabel')) 
  xlabel = deblank(a0.Data.Attributes.MetaData.xlabel);
  xlabel(1:length('# xlabel: '))='';
else xlabel=''; end
if ~isempty(findfield(a0, 'ylabel', 'cache')) 
  ylabel = deblank(a0.Data.Attributes.MetaData.ylabel);
  ylabel(1:length('# ylabel: '))='';
else ylabel=''; end
if ~isempty(findfield(a0, 'xvars', 'cache')) 
  xvars = deblank(a0.Data.Attributes.MetaData.xvars);
  xvars(1:length('# xvars: '))='';
else xvars=''; end

if ~isempty(xvars)
  xvars_i = find(cellfun('isempty', strfind(cnames,xvars)) == 0);
  if ~isempty(xvars_i)
    if length(xvars_i) > 1
      cnames=cnames(xvars_i(end):end);
      xvars_i=xvars_i(1);
    end
    setalias(a0,'x',['this.' getalias(a0,'Signal') '(:,' num2str(xvars_i) ')' ],xvars); % create an alias for xvars
    setalias(a0,xvars,'x',xvars); % create an alias for xvars
    % locate xvars label and column
    xlabel=xvars;
  end

  % Define scanning variable
  setaxis(a0,1,'x');
end


param = load_mcstas_param(a0, 'Param');
a0.Data.Parameters = param;
setalias(a0, 'Parameters', 'Data.Parameters', 'Instrument parameters');

siz = size(a0.Signal);
siz = (siz(2)-1)/2;

a = [];
for j=1:siz
  b = copyobj(a0);
  ylabel=cnames(2*j);
  setalias(b,'Signal', ['this.' getalias(a0,'Signal') '(:,' num2str(2*j) ')'], ylabel);
  if ~isempty(findfield(a0, '_ERR', 'cache')) 
    setalias(b,'Error',['this.' getalias(a0,'Signal') '(:,' num2str(1+2*j) ')']);
  end
  b.Title = [ char(ylabel) ': ' char(b.Title) ];
  b.Label = [ char(ylabel) '(' xvars ')' ];
  a = [a b];
end

% ------------------------------------------------------------------------------
% build-up a parameter structure which holds all parameters from the simulation
function param=load_mcstas_param(a, keyword)
  if nargin == 1, keyword='Param:'; end
  param = [];
  
  par_list = findfield(a, keyword,'case cache');
  if ~isempty(par_list)
    for index=1:numel(par_list)
      try; this=get(a,par_list{index}); catch; this=[]; end
      if isscalar(this) && (isnumeric(this) || ischar(this))
        name         = fliplr(strtok(fliplr(par_list{index}),sprintf('.')));
        param.(name) = this; 
      end
    end
  end

  par_list = findstr(a, keyword, 'case');
  if ischar(par_list), par_list=cellstr(par_list); end
  
  % search strings of the form 'keyword' optional ':', name '=' value
  for index=1:length(par_list)
    line         = par_list{index};
    reversed_line= line(end:-1:1);
    equal_sign_id= find(reversed_line == '=');
    if isempty(equal_sign_id), continue; end % incomplete line: skip it
    name         = fliplr(strtok(reversed_line((equal_sign_id+1):end),sprintf(' \n\t\r\f;#')));
    if isempty(name)
      column_sign_id = findstr(line, keyword);
      name = strtok(line((column_sign_id+length(keyword)+1):end));
    end
    
    if isfield(param, name) % was set with a Previous 'Param name=value' token
      continue;
    elseif isfield(a.Data, name) && any(getfield(a.Data, name))
      value = getfield(a.Data, name);
    else
      value = strtok(fliplr(reversed_line(1:(equal_sign_id-1))),sprintf(' \n\t\r\f;#'));
      if ~isempty(str2num(value)), value = str2num(value); end
    end
    if ~isempty(value) && ~isempty(name) && ischar(name)
      param = setfield(param, name, value);
    end
  end

