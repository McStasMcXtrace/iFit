function obj = copy_prop(obj0, m)
  % copy the properties from 'in' in the instantiated object obj0 (single)
  % handle arrays by copying the obj0
  obj = [];
  for index=1:numel(m)
    this = m(index);

    % transfer properties
    % this is a safe way to instantiate a subclass
    warning off MATLAB:structOnObject
    this = struct(this);
    for p = fieldnames(this)'
      obj0.(p{1}) = this.(p{1});
    end
    obj0.class = class(obj0);
    
    if index == 1
      obj = obj0;
    else
      obj = [ obj obj0 ];
    end
  end
