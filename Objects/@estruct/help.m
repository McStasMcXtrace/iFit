function c = help(s)
% HELP display object information as a Dialogue window or a structure
%   HELP(S) displays the object information in a Dialogue window
%
%   H = HELP(S) returns the object information as a structure
%
% Example: s=estruct(1:10); isstruct(help(s))
%
% Version: $Date$ $Version$ $Author$
% See also  estruct/cell, estruct/double, estruct/struct, 
%           estruct/char, estruct/size, estruct/get
%

% this function is used to gather info for the contextual menu in plots.

c={};
for index=1:numel(s)
  if numel(s) > 1, a=s(index); disp(a); 
  else a=s; end
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

  % info about the plot --------------------------------------------------
  T   = a.Name; ttl = title(a);
  if isempty(T), T=ttl; end
  if ~ischar(T) || isempty(T), T=char(T,'short'); end
  if ~isvector(T), T=transpose(T); T=T(:)'; end
  T   = regexprep(T,'\s+',' '); % remove duplicate spaces
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

  % Label
  d = '';
  if ~isempty(a.Label)
    g = cellstr(a.Label); g=deblank(g{1});
    if length(g) > 23, g = [ g(1:20) '...' ]; end           % Label
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
               [ 'Name: "' char(T) '" ' d ], ...
               [ 'Source: ' a.Source ], ...
               [ 'Last command: ' cmd ]};
  tproperties = {};
  % axes and Signal stuff
  properties{end+1} = '[Rank]         [Value] [Description]';
  % code from estruct.disp -----------------------------------------------------
  myisvector = @(c)length(c) == numel(c);
  for index=0:length(a.Axes)
    [v, l] = getaxis(a, num2str(index,2));
    if ~ischar(v)
      if numel(v) > 5, v=v(1:5); end
      v=mat2str(v);
      if length(v) > 12, v = [v(1:12) '...' ]; end
    end
    if length(l) > 20, l = [l(1:18) '...' ]; end
    X      = getaxis(a, index); x=X(:);
    if issparse(x), x=full(x); end
    if length(x) == 1
      minmaxstd = sprintf('[%g]', x);
    elseif myisvector(X)
      minmaxstd = sprintf('length [%i]', numel(x));
    else
      minmaxstd = sprintf('size %s', mat2str(size(X)));
    end
    t = sprintf('%6i %15s  %s %s\n', index, v, l, minmaxstd);
    tproperties{end+1} = t;
    properties{end+1}  = t;
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
  ud.name   = char(a,'short');
  ud.commands = commandhistory(a);
  ud.handle = [];
  ud.tproperties = tproperties;
  ud.mproperties = mproperties;
  
  % now open the dialogue when not used in nargout=help(object)
  if nargout == 0
    h = helpdlg(ud.properties, [ mfilename ': Data ' a.Tag ': ' num2str(ndims(a)) 'D object ' mat2str(size(a)) ]);
    ud.handle = h;
  end
  
  
  if numel(s) == 1, c=ud; else c{index} = ud; end
end

