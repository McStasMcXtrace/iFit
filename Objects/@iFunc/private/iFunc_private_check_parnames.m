function a = iFunc_private_check_parnames(a, nb_pars)
% iFunc_private_check_parnames: check parameter names
%
%  build missing parameter names
%  check unicity of parameter names
%  replace "parname" with the corresponding 'p(n)' expression

if nargin < 2, nb_pars = []; end

% default parameter names
% try to be clever by an analysis of the expression...
amp=0; cen=0; bkg=0; wid=0;
if ~isempty(nb_pars) && (isempty(a.Parameters) || any(cellfun(@isempty, a.Parameters)))
  expr = a.Expression;
  if isa(expr, 'function_handle'), expr = func2str(expr); end
  if iscellstr(expr), expr=char(expr); end
  for index=1:nb_pars
    namp = [ 'p(' num2str(index) ')' ];
    if ~isempty(findstr(expr, [ '*' namp ])) || ~isempty(findstr(expr, [ namp '*' ]))
      if ~amp
        name = 'Amplitude'; amp=1;
      else
        name = [ 'Amplitude_' num2str(index) ];
      end
    elseif ~isempty(findstr(expr, [ '-' namp ]))
      if ~cen
        name = 'Centre'; cen=0;
      else
        name = [ 'Centre_' num2str(index) ];
      end
    elseif ~isempty(findstr(expr, [ '+' namp ])) || ~isempty(findstr(expr, [ namp '+' ]))
      if ~bkg
        name = 'Constant'; bkg=1;
      else
        name = [ 'Constant_' num2str(index) ];
      end
    elseif ~isempty(findstr(expr, [ '/' namp ])) || ~isempty(findstr(expr, [ namp '^2' ]))  
      if ~wid
        name = 'Width'; wid=1;
      else
        name = [ 'Width_' num2str(index) ];
      end
    else
      name = [ a.Tag '_p' num2str(index) ];
    end
    if numel(a.Parameters) < index || isempty(a.Parameters{index})
      a.Parameters{index} = name;
    end
  end
else
  % re-arrange the parameter names as a string representing a cellstr
  for index=1:length(a.Parameters)
    if ~isempty(a.Parameters{index})
      u = a.Parameters{index};
      u(~isstrprop(u,'print'))=' ';
      a.Parameters{index} = u;
    end
  end
end

% check for unicity of names and possibly rename similar ones
[Pars_uniq, i,j] = unique(strtok(a.Parameters)); % length(j)=Pars_uniq, length(i)=Parameters
for index=1:length(Pars_uniq)
  index_same=find(strcmp(Pars_uniq(index), strtok(a.Parameters)));
  if length(index_same) > 1 % more than one parameter with same name
    for k=2:length(index_same)
      [tok,rem] = strtok(a.Parameters{index_same(k)});
      % check if parameter name is already a renamed one with '_<number>'
      j = regexp(tok, '_[0-9]+$', 'match');
      if ~isempty(j)
        j = j{1};                       % the new incremented parameter duplicate
        j = str2num(j(2:end))+length(index_same)-1; 
        tok((end-length(j)):end) = '';  % the root of the name
      else
        j = k;
      end
      a.Parameters{index_same(k)} = [ tok '_' num2str(j) ' ' rem ];
    end
  end
end

% build the list of replacement strings
replace = strcat('p(', cellstr(num2str(transpose(1:length(a.Parameters)))), ')');
% replace Parameter names by their p(n) representation
if ischar(a.Expression) || iscellstr(a.Expression)
  a.Expression = regexprep(a.Expression, strcat('\<"', a.Parameters, '"\>' ), replace);
end
if ~isempty(a.Constraint.eval) && (ischar(a.Constraint.eval) || iscellstr(a.Constraint.eval))
  a.Constraint.eval = regexprep(a.Constraint.eval, strcat('\<"', a.Parameters, '"\>' ), replace);
end
for index=1:numel(a.Constraint.set)
  if ~isempty(a.Constraint.set{index}) && (ischar(a.Constraint.set{index}) || iscellstr(a.Constraint.set{index}))
    a.Constraint.set{index} = regexprep(a.Constraint.set{index}, strcat('\<"', a.Parameters, '"\>' ), replace);
  end
end
  
