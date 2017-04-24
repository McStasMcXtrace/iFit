function c = help(s)
% c = help(s) : display iData information as a Dialogue window
%
%   @iData/help: display information about the data set as a Dialogue
%
% input:  s: object or array (iData) 
% output: c: iData identification (structure)
%
% Version: $Date$
% See also  iData/cell, iData/double, iData/struct, 
%           iData/char, iData/size, iData/get
%

% EF 23/09/07 iData implementation
% this function is used to gather info for the contextual menu in plots.

c={};
for index=1:numel(s)
  if numel(s) == 1, a=s; else a=s(index); end
  % Model stuff ----------------------------------------------------------------
  m = []; mp = []; mv = []; names = []; name = [];
  % get Model,etc... when found in the Dataset

  if isfield(a, 'Model')
    m = get(a, 'Model');
  elseif ~isempty(findfield(a, 'Model'))
    m = get(a, findfield(a, 'Model', 'cache first'));
  end

  if isa(m, 'iFunc') && ~isempty(m)
    % get the parameter values as a struct
    mp    = m.ParameterValues;
    names = m.Parameters; names = names(:);
    name  = m.Name;
  end


  if isfield(a, 'ModelValue')
    mv = get(a, 'ModelValue');
    if ~strcmp(getalias(a,'Signal'), 'ModelValue')
      set(a, 'ModelValue', []);  % avoid recursive loop
    end
  elseif ~isempty(findfield(a, 'ModelValue'))
    mv = get(a, findfield(a, 'ModelValue', 'cache first'));
  end

  if isempty(mp)
    if isfield(a, 'ModelParameters')
      mp = get(a, 'ModelParameters');
    elseif isfield(mv, 'ModelParameters')
      mp = get(mv, 'ModelParameters');
    end
  end
  % make it a structure
  if ~isempty(mp) && numel(names) == numel(mp)
    mp = cell2struct(num2cell(mp(:)),strtok(names(:)));
  end

  % info about the iData plot --------------------------------------------------
  T   = a.Title; if ~ischar(T), T=char(T); end
  if ~isvector(T), T=transpose(T); T=T(:)'; end
  T   = regexprep(T,'\s+',' '); % remove duplicated spaces
  cmd = char(a.Command{end});
  S   = a.Source;
  [pS, fS, eS] = fileparts(S);
  if length(pS) > 13, pS=[ '...' pS(end-10:end) ]; end
  if length(fS) > 13, fS=[ '...' fS(end-10:end) ]; end
  if ~isempty(pS), S = [ pS filesep ];
  else             S = '';
  end
  S = [ S fS ];
  if ~isempty(eS), S = [ S '.' eS ]; end
  if length(cmd) > 23, cmd = [ cmd(1:20) '...' ]; end

  % DisplayName and Label
  d = '';
  if ~isempty(a.Label) && ~isempty(a.DisplayName)
    if strcmp(a.Label, a.DisplayName)
        if ~isempty(title(a)), a.DisplayName=title(a);
        else a.DisplayName=fS; end
    end
    g = cellstr(a.Label); g=deblank(g{1});
    if length(g) > 13, g = [ g(1:10) ]; end                 % Label/DisplayName
    d = [ d sprintf('%s', g) ];
    g = cellstr(a.DisplayName); g=deblank(g{1});
    if length(g) > 13, g = [ g(1:10) '...' ]; end           % 
    d = [ d sprintf('/%s', g) ];
  elseif ~isempty(a.Label)
    g = cellstr(a.Label); g=deblank(g{1});
    if length(g) > 23, g = [ g(1:20) '...' ]; end           % Label
    d = [ d sprintf('%s', g) ];
  elseif ~isempty(a.DisplayName)
    g = cellstr(a.DisplayName); g=deblank(g{1});
    if length(g) > 23, g = [ g(1:20) '...' ]; end           % DisplayName
    d = [ d sprintf('%s', g) ];
  end
  T0 = T; % original title, full.

  if length(T) > 23, T=[ T(1:20) '...' ]; end
  if length(S)+length(d) < 30,
    d = [ d ' ' T ];
  end
  
  % ----------------------------------------------------------------------------
  % make up title string and Properties dialog content
  properties={ [ 'Data ' a.Tag ': ' num2str(ndims(a)) 'D object ' mat2str(size(a)) ], ...
               [ 'Title: "' char(T) '" ' d ], ...
               [ 'Source: ' a.Source ], ...
               [ 'Last command: ' cmd ]};
  tproperties = {};
  % axes and Signal stuff
  properties{end+1} = '[Rank]         [Value] [Description]';
  myisvector = @(c)length(c) == numel(c);
  for index=0:min([ ndims(a) length(getaxis(a)) ])
    [v, l] = getaxis(a, num2str(index));
    if length(l) > 20, l = [l(1:18) '...' ]; end 
    x      = getaxis(a, index);
    m      = get(a, 'Monitor');
    if length(x) == 1
      minmaxstd = sprintf('[%g]', full(x));
    elseif myisvector(x)
      minmaxstd = sprintf('[%g:%g] length [%i]', full(min(x)), full(max(x)),length(x));
    else
      x=x(:);
      minmaxstd = sprintf('[%g:%g] size [%s]', full(min(x)), full(max(x)),num2str(size(x)));
    end
    if index==0
      if not(all(m==1 | m==0))
        minmaxstd=[ minmaxstd sprintf(' (per monitor=%g)', mean(m(:))) ];
      end
      minmaxstd=[ minmaxstd sprintf(' sum=%g', full(sum(iData_private_cleannaninf(x)))) ];
    end
    if prod(size(a)) < 1e4
      try
        [s, f] = std(a, -index);
        minmaxstd=[ minmaxstd sprintf(' <%g +/- %g>', f,s) ];
      end
    end
    if isnumeric(v), v=''; end
    t = sprintf('%6i %15s  %s %s', index, v, l, minmaxstd);
    tproperties{end+1} = t;
    properties{end+1}  = t;
    clear x m
  end

  % model parameters
  if ~isempty(mp)
    mproperties = { ['Model parameters: ' name ] };
    if isstruct(mp)
      for f=fieldnames(mp)'
        mproperties{end+1} = sprintf('* %s = %g', f{1}, mp.(f{1}));
      end
    elseif isnumeric(mp)
      mproperties{end+1} = mat2str(mp);
    end
  else mproperties = {};
  end
  properties = { properties{:} ...
     ' ' ...
     mproperties{:} };

  % attach contexual menu to plot with UserData storage
  ud.properties=properties;
  ud.xlabel = xlabel(a);
  ud.ylabel = ylabel(a);
  ud.zlabel = zlabel(a);
  ud.title  = T;
  ud.name   = char(a);
  ud.commands = commandhistory(a);
  ud.handle = [];
  ud.tproperties = tproperties;
  ud.mproperties = mproperties;

  % ud = iData_plot_contextmenu(a, [], xlabel(a), ylabel(a), zlabel(a),  T, S, d, cmd, mp, name);
  
  % now open the dialogue when not used in nargout=help(object)
  if nargout == 0
    h = helpdlg(ud.properties, [ mfilename ': Data ' a.Tag ': ' num2str(ndims(a)) 'D object ' mat2str(size(a)) ]);
    ud.handle = h;
  end
  
  
  if numel(s) == 1, c=ud; else c{index} = ud; end
end

