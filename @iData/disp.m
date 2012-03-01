function disp(s_in, name)
% disp(s) : display iData object (details)
%
%   @iData/disp function to display iData object details
%
% input:  s: object or array (iData) 
% ex:     'disp(iData)'
%
% Version: $Revision: 1.26 $
% See also iData, iData/display, iData/get

% EF 27/07/00 creation
% EF 23/09/07 iData implementation

if nargin == 2 && ~isempty(name)
  iname = name;
elseif ~isempty(inputname(1))
  iname = inputname(1);
else
  iname = 'ans';
end

% removes warnings during disp
iData_private_warning('enter',mfilename);
 
if numel(s_in) > 1
  eval([ iname ' = s_in;' ])
  eval([ 'display(' iname ');' ]); % makes sure the variable name is sent to 'display'.
else
  fprintf(1,'%s = iData %iD object of size [%s]:\n',iname, ndims(s_in), num2str(size(s_in)));
  m = get(s_in, 'Monitor'); m=m(:);
  s=struct(s_in);
  s=rmfield(s,'Alias');
  if isnumeric(s.Date), s.Date=datestr(s.Date); end
  if isnumeric(s.ModificationDate), s.ModificationDate=datestr(s.ModificationDate); end
  disp(s)
  % display the Aliases
  disp('Object aliases:');
  disp('         [Name]                           [Value]   [Description]');
  for index=1:length(s_in.Alias.Names)
    v = s_in.Alias.Values{index};
    if isempty(v)
      if strcmp(s_in.Alias.Names{index}, 'Error'), v='[] i.e. sqrt(Signal)';
      elseif strcmp(s_in.Alias.Names{index}, 'Monitor'), v='1'; end
    end 
    label = s_in.Alias.Labels{index};
    if length(label) > 30, label = [label(1:28) '...' ]; end 
    if isnumeric(v) | islogical(v), 
      if length(size(v)) > 2, v=v(:); end
      if numel(v) > 20, v=v(1:18); end
      v = mat2str(v,2); 
      if length(v) > 32, v = [v(1:29) '...' ]; end 
    elseif ischar(v)
      this = s_in;
      try
        vv = eval(v);
        if length(vv) > 20, vv = vv(1:18); end 
        vv = mat2str(vv);
        if length(vv) > 20, vv = [vv(1:18) '...' ]; end 
        label = [ label ' ''' vv '''' ];
      catch
        try
          vv = get(s_in, v);
          if isnumeric(vv), 
            if length(size(vv)) > 2, vv=vv(:); end
            if numel(vv) > 10, vv=vv(1:10); end
            vv = mat2str(vv,2); 
          elseif ischar(vv)
            vv = vv(:)';
          end
          if length(vv) > 20, vv = [vv(1:18) '...' ]; end
          label = [ label ' ''' vv '''' ];
        end
      end
      if length(v) > 32, v = [ '...' v((end-28):end) ]; end
    end
    fprintf(1,'%15s  %32s   %s', s_in.Alias.Names{index}, strtrim(v), strtrim(char(label)));
    if strcmp(s_in.Alias.Names{index}, 'Signal') & length(s_in.Alias.Axis) == 0
      x      = getaxis(s_in, 0);  x=x(:);
      if length(x) == 1
        fprintf(1,' [%g]', x);
      else
        fprintf(1,' [%g:%g]', min(x), max(x));
      end
      if ~(all(m(:)==1) | all(m(:)==0))
        fprintf(1,' (per monitor)');
      end
    end
    fprintf(1, '\n');
  end
  
  % display the Signal and Axes
  disp('Object axes:');
  disp('[Rank]         [Value]  [Description]');
  for index=0:length(s_in.Alias.Axis)
    [v, l] = getaxis(s_in, num2str(index,2));
    if ~ischar(v)
      if numel(v) > 5, v=v(1:5); end
      v=mat2str(v);
      if length(v) > 12, v = [v(1:12) '...' ]; end 
    end
    if length(l) > 20, l = [l(1:18) '...' ]; end 
    X      = getaxis(s_in, index); x=X(:);
    if length(x) == 1
      minmaxstd = sprintf('[%g]', full(x));
    elseif isvector(X)
      minmaxstd = sprintf('[%g:%g] length [%i]', full(min(x)), full(max(x)),length(x));
    else
      minmaxstd = sprintf('[%g:%g] size [%s]', full(min(x)), full(max(x)),num2str(size(X)));
    end
    if index==0
      if not(all(m==1 | m==0))
        minmaxstd=[ minmaxstd sprintf(' (per monitor=%g)', mean(m(:))) ];
      end
      minmaxstd=[ minmaxstd sprintf(' sum=%g', sum(iData_private_cleannaninf(x))) ];
    end
    if prod(size(s_in)) < 1e4
      try
        [s, f] = std(s_in, index);
        minmaxstd=[ minmaxstd sprintf(' <%g +/- %g>', f,s) ];
      end
    end
    fprintf(1,'%6i %15s  %s %s\n', index, v, l, minmaxstd);
  end
  %iData(s_in);
end

% reset warnings during disp
iData_private_warning('exit',mfilename);

