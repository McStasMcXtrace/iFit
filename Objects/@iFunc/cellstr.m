function [ret, header] = cellstr(s)
% ret = cellstr(s) : convert iFunc into a cell of strings
%
%   @iFunc/cellstr: function to convert iFunc objects into a cell of
%       strings e.g. for further evaluation. 
%   Returns the iFunc expression to evaluate
%
% input:  s: object or array (iFunc) 
% output: ret: iFunc identification (cellstr)
%
% Version: $Revision: 1.1 $
% See also  iFunc/struct, iFunc/char
%


ret=[];
if numel(s) > 1
  ret = {};
  for index=1:numel(s)
    ret{index} = char(s(index));
  end
  return
end
  
% single object to char
  if isempty(s), ret='[]'; return; end

  ax = 'x,y,z,t,u,'; ax = ax(1:(s.Dimension*2));
  
  if isempty(s.Name) || strcmp(s.Name, char(s.Expression)), n  = s.Tag; else n=s.Name; end
  if strcmp(s.Expression, s.Description),             d = '';     else d = s.Description; end
  
  NL = sprintf('\n');
  % now we build up the header
  ret = { sprintf('%% signal=%s(p,%s ...) iFunc object %s', n, ax, s.Tag), '%' };
  if ~isempty(d) 
    d = textwrap(cellstr(d),80); 
    for index=1:length(d)
      ret{end+1} = sprintf('%% %s', d{index});
    end
    ret{end+1} = '%';
  end
  
  ret{end+1} = '% Expression:';
  if ~isempty((s.Constraint)) 
    e = textwrap(cellstr(s.Constraint),80);
    if length(e) > 3, e=e(1:3); e{end+1} = '...'; end
    for index=1:length(e)
      ret{end+1} = sprintf('%%   %s', e{index});
    end
  end
  e = textwrap(cellstr(char(s.Expression)),80);
  if length(e) > 3, e=e(1:3); e{end+1} = '...'; end
  for index=1:length(e)
    ret{end+1} = sprintf('%%   %s', e{index});
  end
  ret{end+1} = '%';
  ret{end+1} = '% input:';
  
  ret{end+1} = [ '%   p: Parameters (' num2str(length(s.Parameters)) ' values, degrees of freedom)' ];
  p = cellstr(s.Parameters);
  for index=1:length(p)
    ret{end+1} = sprintf('%%      p(%2d)=%s', index, p{index});
  end
  ret{end+1} = [ '%   ' ax(1:(end-1)) ': model axes (' num2str(s.Dimension) ' vector/matrix, dimensionality)' ];
  ret{end+1} = '%   ...: additional arguments to the model function';
  ret{end+1} = '% output:';
  ret{end+1} = '%   signal: function value or information';
  ret{end+1} = '';
  header = char(ret);

  % now write the core of the model (for evaluation)
  if ~isempty(s.Constraint) 
    if isa(s.Constraint ,'function_handle')
      ret{end+1} = sprintf('p2 = feval(%s, p, %s); p(~isnan(p2))=p2(~isnan(p2));', fun2str(s.Constraint), ax(1:(end-1)));
    else
      ret{end+1} = '% The Constraint';
      e = cellstr(s.Constraint);
      for index=1:length(e)
        this = strtrim(e{index});
        if this(end) == ';'
          ret{end+1} = sprintf('%s\n', e{index});
        else
          ret{end+1} = sprintf('%s;\n', e{index});
        end
      end
    end
  end
  
  % this one has to return a 'signal' value in the last line
  if isa(s.Expression ,'function_handle')
    ret{end+1} = sprintf('signal = feval(%s, p, %s);', func2str(s.Expression), ax(1:(end-1)));
  else
    ret{end+1} = [ '% The Expression, computing signal from ''p'' and axes ' ax(1:(end-1)) ];
    e = s.Expression; 
    if ischar(e) % split char into lines
      e = textscan(e,'%s','Delimiter',sprintf('\n\r\f'),'MultipleDelimsAsOne',1); e=e{1};
    end
    has_signal = 0;
    for index=1:length(e)
      d = strtrim(e{index});
      if d(end) ~= ';', d = [ d '; ' ]; end
      if ~isempty(regexp(d, '\<signal\>\s*=')), has_signal = 1; end
      if index == length(e) && ~has_signal
        ret{end+1} = sprintf('signal = %s', d);
        has_signal = 1;
      else
        ret{end+1} = sprintf('%s', d);
      end
    end
  end

  % return value
  ret = reshape(ret,numel(ret),1); % as rows

  if nargout == 0 && ~isempty(inputname(1))
    s.Eval = ret;
    assignin('caller',inputname(1),s);
  end

