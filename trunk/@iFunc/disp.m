function disp(s_in, name)
% disp(s) : display iFunc object (details)
%
%   @iFunc/disp function to display iFunc model details
%
% input:  s: object or array (iFunc) 
% ex:     'disp(iFunc)'
%
% Version: $Revision: 1.2 $
% See also iFunc, iFunc/display, iFunc/get

if nargin == 2 && ~isempty(name)
  iname = name;
elseif ~isempty(inputname(1))
  iname = inputname(1);
else
  iname = 'ans';
end

if numel(s_in) > 1
  eval([ iname ' = s_in;' ])
  eval([ 'display(' iname ');' ]); % makes sure the variable name is sent to 'display'.
else
  if isdeployed || ~usejava('jvm'), id='iFunc';
  else           id='<a href="matlab:helpwin iFunc">iFunc</a>';
  end
  fprintf(1,'%s = %s %iD model:\n',iname, id, s_in.Dimension);
  % clean up redundant/empty fields
  s = struct(s_in);

  u = char(s.Constraint); u=strtrim(u); u(~isstrprop(u,'print'))=''; if ~isvector(u), u=u'; end
  if length(u) > 70, u = [ u(1:67) '...' ]; end
  if ~isempty(u)
    fprintf(1, '         Constraint: %s\n', u); 
  end
  u = char(s.Expression); u=strtrim(u); u(~isstrprop(u,'print'))=''; if ~isvector(u), u=u'; end
  if length(u) > 70, u = [ u(1:67) '...' ]; end
  fprintf(1, '         Expression: %s\n', u); 
  if strcmp(s.Name, s.Description)       || isempty(s.Name),        s =rmfield(s, 'Name'); end
  if strcmp(s.Expression, s.Description) || isempty(s.Description), s =rmfield(s, 'Description'); 
  else
    u = s.Description; u=strtrim(u); u(~isstrprop(u,'print'))=''; if ~isvector(u), u=u'; end
    fprintf(1, '        Description: %s\n', s.Description); 
    s =rmfield(s, 'Description');
  end
  if isempty(s.Guess),                    s =rmfield(s, 'Guess'); end
  if isempty(s.Constraint),               s =rmfield(s, 'Constraint'); end
  s=rmfield(s, 'Eval');
  s=rmfield(s, 'Expression');
  if isnumeric(s.Date), s.Date=datestr(s.Date); end
  % object Properties displayed as a structure
  disp(s)
  % now display parameters in compact form
  if ~isempty(s.Parameters)
    disp('Parameters:')
    for p=1:length(s.Parameters)
      [name, R] = strtok(s.Parameters{p}); % make sure we only get the first word (not following comments)
      R = strtrim(R);
      fprintf(1,'  p(%3d)=%20s', p, name);
      val  = [];
      if ~isempty(s.ParameterValues)
        try
          val = s.ParameterValues(p);
        end
      end
      if ~isempty(val), fprintf(1, '=%g', val); end
      if ~isempty(R),   fprintf(1, '  %% entered as: %s', R); end
      fprintf(1, '\n');
    end

  end
  
end


