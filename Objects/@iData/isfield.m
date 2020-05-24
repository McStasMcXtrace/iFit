function tf = isfield(self, f)
  % ISFIELD True if field is in object.
  %   ISFIELD(S,FIELD) returns true if the string FIELD is the name of a
  %   field in the object S. FIELD must be a single char.
  %
  %   NOTE: TF is false when FIELD or FIELDNAMES are empty.
  %
  %   Example: s = iData('one',1,'two',2); isfield(s,'one')
  %
  % See also getfield, setfield, fieldnames, orderfields, rmfield,
  %   isstruct, iData.
    tf=false;
    if nargin < 2, f=''; end
    % handle case for multiple field names
    if numel(self) == 1 && iscellstr(f)
      tf = logical(zeros(size(f)));
      for index=1:numel(f)
        tf(index) = isfield(self, f{index});
      end
      return
    end
    
    if isempty(f) || (~ischar(f) && ~iscellstr(f)), return; end
    if numel(self) == 1
      if any(f == '.') % compound field: we search with findfield
        tf = findfield(self, f, 'isfield');
      else
        tf = any(strcmp(f, fieldnames(self)));
        % tf = isfield(struct(self), f);
        % tf = ~isempty(findprop(self, f));
      end
    else
      tf = cell(size(self)); is_scalar = true;
      for index=1:numel(self);
        t = isfield(self(index), f);
        if ~isscalar(t), is_scalar=false; end
        tf{index} = t;
      end
      if is_scalar, tf = cell2mat(tf); end
    end
  end
