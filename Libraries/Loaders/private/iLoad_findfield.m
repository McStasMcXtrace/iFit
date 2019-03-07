function match = findfield(s, field)
% match=findfield(s, field, option) : look for numerical fields in a structure
%
% input:  s: structure
%         field: field name to search, or '' (char).
%         option: 'exact' 'case' or '' (char)
% output: match: names of structure fields (cellstr)
% ex:     findfield(s) or findfield(s,'Title')

  if numel(s) > 1
    match = cell(1, numel(s));
    for index=1:length(s)
      match{index}=findfield(s(index), field);
    end
    return
  end

  if iscell(s), s=s{1}; end
  match = struct_getfields(struct(s), ''); % return the numerical fields

  if ~isempty(field)
    field = lower(field);
    matchs= lower(match);

    if iscellstr(field)
      index = [];
      for findex=1:length(field)
        tmp = strfind(matchs, field{findex});
        if iscell(tmp), tmp = find(cellfun('isempty', tmp) == 0); end
        index= [ index ; tmp ];
      end
      index = unique(index);
    else
      index = strfind(matchs, field);
    end

    if ~isempty(index) && iscell(index), index = find(cellfun('isempty', index) == 0); end
    if isempty(index)
      match=[];
    else
      match = match(index);
    end
  end

end % findfield

% ============================================================================
% private function struct_getfields, returns field, class, numel 
function f = struct_getfields(structure, parent)

  f=[];
  if ~isstruct(structure), return; end
  if numel(structure) > 1
    structure=structure(:);
    for index=1:length(structure)
      sf = struct_getfields(structure(index), [ parent '(' num2str(index) ')' ]);
      f = [f(:) ; sf(:)];
    end
    return
  end

  % get content and type of structure fields
  c = struct2cell(structure);
  f = fieldnames(structure);
  try
    t = cellfun(@class, c, 'UniformOutput', 0);
  catch
    t=cell(1,length(c));
    index = cellfun('isclass', c, 'double'); t(find(index)) = {'double'};
    index = cellfun('isclass', c, 'single'); t(find(index)) = {'single'};
    index = cellfun('isclass', c, 'logical');t(find(index)) = {'logical'};
    index = cellfun('isclass', c, 'struct'); t(find(index)) = {'struct'};
    index = cellfun('isclass', c, 'uint8');  t(find(index)) = {'uint8'};
    index = cellfun('isclass', c, 'uint16'); t(find(index)) = {'uint16'};
    index = cellfun('isclass', c, 'uint32'); t(find(index)) = {'uint32'};
    index = cellfun('isclass', c, 'uint64'); t(find(index)) = {'uint64'};
    index = cellfun('isclass', c, 'int8');   t(find(index)) = {'int8'};
    index = cellfun('isclass', c, 'int16');  t(find(index)) = {'int16'};
    index = cellfun('isclass', c, 'int32');  t(find(index)) = {'int32'};
    index = cellfun('isclass', c, 'int64');  t(find(index)) = {'int64'};
  end

  toremove=[];
  % only retain numerics
  for index=1:length(c)
    if ~any(strncmp(t{index},{'dou','sin','int','uin','str','log'}, 3))
      toremove(end+1)=index;
    end
  end
  c(toremove)=[];
  f(toremove)=[];

  if ~isempty(parent), f = strcat([ parent '.' ], f); end

  % find sub-structures and make a recursive call for each of them
  for index=transpose(find(cellfun('isclass', c, 'struct')))
    try
    sf = struct_getfields(c{index}, f{index});
    f = [f(:) ; sf(:)];
    end
  end

end % struct_getfields
