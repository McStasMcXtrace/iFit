function ret=private_id
% PRIVATE_ID generate a new unique ID

  persistent id
  
  % handle Tag number
  if isempty(id) id=0; end
  if id > 1e6,   id=0; end % use clock
  if id <=0,
    id = clock;
    id = fix(id(6)*1e4);
  else
    id=id+1;
  end

  ret = id;
