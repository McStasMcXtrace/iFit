function tf = isfield(self, f)
  % ISFIELD True if field is in object.
  %   ISFIELD(S,FIELD) returns true if the string FIELD is the name of a
  %   field in the object S. FIELD must be a single char.
  %
  %   NOTE: TF is false when FIELD or FIELDNAMES are empty.
  %
  %   Example: s = estruct('one',1,'two',2); isfield(s,'one')
  %
  % See also getfield, setfield, fieldnames, orderfields, rmfield,
  %   isstruct, estruct.
    tf=false;
    if nargin < 2, f=''; end
    if isempty(f) || (~ischar(f) && ~iscellstr(f)), return; end
    if numel(self) == 1
      if any(f == '.') % compound field: we search with findfield
        tf = findfield(self, f, 'isfield');
      else
        tf = isfield(struct(self), f); % faster
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