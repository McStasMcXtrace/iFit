function b = subsref(a,S)
% b = subsref(a,s) : Process indexed references
%
%   @Process/subsref: function returns subset indexed references
%     such as a(1:2) or a.field.
%             a{n} returns the guessed axis of rank 'n' (similar to iData).
%
% Version: $Date$
% See also Process, Process/subsasgn

% This implementation is very general, except for a few lines

persistent fields

if isempty(fields), fields=fieldnames(Process); end

b = a;  % will be refined during the index level loop

if isempty(S)
  return
end

for i = 1:length(S)     % can handle multiple index levels
  s = S(i);
  switch s.type
  case '()' % ======================================================== array
    if numel(b) > 1           % Process array: b(index)
      b = b(s.subs{:});
    end
  case '.'  % ======================================================== structure
    % protect some fields     % Process Property
    fieldname = s.subs;
    if length(fieldname) > 1 && iscell(fieldname)
      fieldname = fieldname{1};
    end
    if isa(b, 'Process'), f=fields; else f=fieldnames(b); end
    index = find(strcmpi(fieldname, f));
    if ~isempty(index) % structure/class def fields: b.field
      b = b.(f{index});
      if isnumeric(b) && (strcmpi(fieldname, 'creationDate') || strcmpi(fieldname, 'terminationDate'))
        b = datestr(b);
      end
    elseif ismethod(b, fieldname)
      if i == length(S)
        if nargout(fieldname) ==0
          builtin('feval',fieldname, b);
          c=[];
        else
          c = builtin('feval',fieldname, b);
        end
      else
        c = builtin('feval',fieldname, b, S(i+1).subs{:}); i=i+1;
      end
      if isa(c, 'Process'), b = c; end
    else
      error([ mfilename ': can not get Process object Property ''' fieldname ''' in Process model ' char(b.process) '.' ]);
    end
  end   % switch s.type
end % for s index level 
