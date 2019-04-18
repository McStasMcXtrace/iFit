function new = copyobj(self, org)
  % COPYOBJ makes a deep copy of initial object.
  %   COPYOBJ(A) copy the properties from A into a new object.
  %   The notation +A is equivalent.
  %
  %   new = COPYOBJ(self) creates a deep copy of the initial object 'self'
  %
  %   new = COPYOBJ(self, content) creates an empty object and fill it
  %   with 'content' (struct/object)
  %
  % Example: s = estruct(1:10); s1=copyobj(s); axescheck(s1); isequal(s.Signal, s1.Signal)
  % Version: $Date$ $Version$ $Author$
  % see also estruct

  % handle arrays by copying the new0
  if nargin == 1,  org=''; end

  if ~isempty(org) && ~isstruct(org) && ~isobject(org)
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
  if numel(org) > 1 && nargin > 1
    new = [];
    for index=1:numel(org)
      new = [ new copyobj(self, org(index)) ];
    end
    return
  end

  % single deep copy
  if ~isa(self, 'handle') && isa(org, 'estruct') && isempty(self)
    % assigment is OK for non handle objects
    new = org;
  else
    % handle object: transfer properties: this is a safe way to instantiate a subclass
    new = feval(class(self)); % new empty object of same class
    wrn_flag = true;
    if isempty(org)
      flag_nargin=1;
      org = self;
    else
      flag_nargin=nargin;
    end
    for p = fieldnames(org)'
      % skip Tag that must be unique
      if strcmp(p{1}, 'Tag'), continue; end
      try
        if ~isfield(new, p{1})
          new.addprop(p{1});
        end
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
  if flag_nargin==1
    history(new, mfilename, self);
  else
    history(new, mfilename, self, org);
  end