function disp(s_in)
% disp(s) : display iData object (details)
%
%   @iData/disp function to display iData object details
%
% input:  s: object or array (iData) 
% ex:     'disp(iData)'
%
% Version: $Revision: 1.7 $
% See also iData, iData/display, iData/get

% EF 27/07/00 creation
% EF 23/09/07 iData implementation

if ~isempty(inputname(1))
  iname = inputname(1);
else
  iname = 'ans';
end
 
if length(s_in) > 1
  eval([ iname ' = s_in;' ])
  eval([ 'display(' iname ');' ]); % makes sure the variable name is sent to 'display'.
else
  fprintf(1,'%s = iData %iD object of size [%s]:\n',iname, ndims(s_in), num2str(size(s_in)));
  s=struct(s_in);
  s=rmfield(s,'Alias');
  disp(s)
  disp('Object aliases:');
  disp('          [Name]                           [Value] [Description]');
  for index=1:length(s_in.Alias.Names)
    v = s_in.Alias.Values{index};
    if isempty(v)
      if strcmp(s_in.Alias.Names{index}, 'Error'), v='sqrt(Signal)';
      elseif strcmp(s_in.Alias.Names{index}, 'Monitor'), v='1'; end
    end 
    label = s_in.Alias.Labels{index};
    if isnumeric(v) | islogical(v), 
      if length(size(v)) > 2, v=v(:); end
      if numel(v) > 10, v=v(1:10); end
      v = mat2str(v); 
      if length(v) > 15, v = [v(1:12) '...' ]; end 
    elseif ischar(v)
      this = s_in;
      try
        vv = mat2str(eval(v));
        if length(vv) > 20, vv = [vv(1:18) '...' ]; end 
        label = [ label ' ''' vv '''' ];
      catch
      end
    end
    fprintf(1,'%15s  %32s   %s', s_in.Alias.Names{index}, v, label);
    if strcmp(s_in.Alias.Names{index}, 'Signal') & length(s_in.Alias.Axis) == 0
      x      = get(s_in, 'Signal');  x=x(:);
      m      = get(s_in, 'Monitor'); m=m(:);
      if ~(all(m==1) | all(m==0))
        x=x./m;
      end
      if length(x) == 1
        fprintf(1,' [%g]', x);
      else
        fprintf(1,' [%g:%g]', min(x), max(x));
      end
      if ~(all(m==1) | all(m==0))
        fprintf(1,' (per monitor)\n');
      else
        fprintf(1,'\n');
      end
    end
    fprintf(1, '\n');
  end
  if length(s_in.Alias.Axis)
    disp('Object axes:');
    disp('[Rank]         [Value] [Description]');
    for index=0:length(s_in.Alias.Axis)
      [v, l] = getaxis(s_in, num2str(index));
      x      = getaxis(s_in, index); x=x(:);
      m      = get(s_in, 'Monitor'); m=m(:);
      if ~(all(m==1) | all(m==0)) & index==0
        x=x./m;
      end
      if length(x) == 1
        fprintf(1,'%6i %15s  %s [%g]', index, v, l, x);
      else
        fprintf(1,'%6i %15s  %s [%g:%g]', index, v, l, min(x), max(x));
      end
      if index==0 & not(all(m==1) | all(m==0))
        fprintf(1,' (per monitor)\n');
      else
        fprintf(1,'\n');
      end
    end
  end
  iData(s_in);
end
