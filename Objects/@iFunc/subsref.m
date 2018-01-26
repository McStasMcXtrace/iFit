function b = subsref(a,S)
% b = subsref(a,s) : iFunc indexed references
%
%   @iFunc/subsref: function returns subset indexed references
%     such as a(1:2) or a.field.
%             a{n} returns the guessed axis of rank 'n' (similar to iData).
%
%   a.p     returns the parameter values, same as a.ParameterValues
%   a.Par   returns the parameter "Par" value
%
% Version: $Date$
% See also iFunc, iFunc/subsasgn

% This implementation is very general, except for a few lines

persistent fields

if isempty(fields), fields=fieldnames(iFunc); end

b = a;  % will be refined during the index level loop

if isempty(S)
  return
end

for i = 1:length(S)     % can handle multiple index levels
  s = S(i);
  switch s.type
  case '()' % ======================================================== array
    if numel(b) > 1           % iFunc array: b(index)
      b = b(s.subs{:});
    elseif isa(b,'iFunc')     % syntax iFunc(p, axes{:}, varargin) -> evaluate
      [b, model, ax, name] = feval(b, s.subs{:});
      % update object
      if nargout == 0 && ~isempty(inputname(1))
        assignin('caller',inputname(1),model);
      end
    else
      b=subsref(b, s); return;
    end
  case '{}' % ======================================================== axes (guessed)
    [~,~,ax] = feval(b);
    b = ax{s.subs{:}};
  case '.'  % ======================================================== structure
    % protect some fields     % iFunc Property
    fieldname = strtok(s.subs);
    if length(fieldname) > 1 && iscell(fieldname)
      fieldname = fieldname{1};
    end
    if isa(b, 'iFunc'), f=fields; else f=fieldnames(b); end
    index = find(strcmpi(fieldname, f));
    if ~isempty(index) % structure/class def fields: b.field
      b = b.(f{index});
      if all(isnumeric(b)) && all(strcmpi(fieldname, 'Date'))
        b = datestr(b);
      end
    elseif (isa(b,'iFunc') || isfield(b, 'Parameters')) && iscellstr(b.Parameters) && any(strcmp(fieldname, strtok(b.Parameters))) % b.<parameter name>
      index=find(strcmp(fieldname, strtok(b.Parameters)));
      if index <= length(b.ParameterValues) % last parameter value used
        b = b.ParameterValues;
        b = b(index);
      else
        b = [];
      end
    elseif strcmpi(fieldname, 'p')               % b.p
      if ~isempty(b.ParameterValues)
        b = b.ParameterValues;
      else
        b = [];
      end
    elseif strcmpi(fieldname, 'Title')           % b.Name
      b = b.Name;
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
      if isa(c, 'iFunc'), b = c; end
    elseif ~isempty(superclasses([ 'iFunc' fieldname ]))
      b = feval([ 'iFunc' fieldname ], b);  % convert to iFunc_<flavour>
    else
      if isa(b, 'iFunc')
        error([ mfilename ': can not get iFunc object Property ''' fieldname ''' in iFunc model ' b.Tag '.' ]);
      else
        error([ mfilename ': can not get Property ''' fieldname ''' in ' class(b) '.' ]);
      end
    end
  end   % switch s.type
end % for s index level 
