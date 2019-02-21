function new = copyobj(self, org)
  % COPYOBJ makes a deep copy of initial object
  %
  %   new = copyobj(self)
  %     creates a deep copy of the initial object 'self'
  %   new = copyobj(self, content)
  %     creates an empty object and fill it with 'content' (struct/object)
  %
  % copy the properties from 'self' into the instantiated object 'new'.
  
  % handle arrays by copying the new0
  if nargin == 1,  org=''; end
  if isempty(org), org=self; end
  if ~isstruct(org) || ~isobject(org)
    new = self; % this may not be a deep copy
    return
  end
  
  % handle array input
  if numel(self) > 1
    new = [];
    for index=1:numel(self)
      new = [ new copyobj(self(index), org) ];
    end
    return
  end
  if numel(org) > 1
    new = [];
    for index=1:numel(org)
      new = [ new copyobj(self, org(index)) ];
    end
    return
  end

  % single deep copy
  if isstruct(self) && ~isa(self, 'handle')
    % assigment is OK for non handle objects
    new = org;  
  else 
    % handle object: transfer properties: this is a safe way to instantiate a subclass
    new = feval(class(self)); % new empty object of same class
    wrn_flag = true;
    for p = properties(org)'
      try
        new.(p{1}) = org.(p{1});  % may fail when copying from enriched object
      catch ME
        if wrn_flag
          disp(getReport(ME))
          disp([ mfilename ': can not copy property ' p{1} ...
                    ' from class ' class(org) ' into ' class(self) ]);
          wrn_flag = false;
        end
      end
    end
  end
  
