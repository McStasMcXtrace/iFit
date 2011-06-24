function b = subsref(a,S)
% b = subsref(a,s) : iData indexed references
%
%   @iData/subsref: function returns subset indexed references
%     such as a(1:2) or a.field.
%   The special syntax a{0} where a is a single iData returns the 
%     Signal/Monitor, and a{n} returns the axis of rank n.
%
% Version: $Revision: 1.14 $
% See also iData, iData/subsasgn

% This implementation is very general, except for a few lines
% EF 27/07/00 creation
% EF 23/09/07 iData impementation

b = a;  % will be refined during the index level loop

if isempty(S)
  return
end

for i = 1:length(S)     % can handle multiple index levels
  s = S(i);
  switch s.type
  case '()'             % extract Data using indexes
    if length(b(:)) > 1   % iData array
      b = b(s.subs{:});
    else                  % single iData
      try % disable some warnings
        warn.seta = warning('off','iData:setaxis');
        warn.geta = warning('off','iData:getaxis');
        warn.get  = warning('off','iData:get');
      catch
        warn = warning('off');
      end
      % this is where specific class structure is taken into account
      if any(cellfun('isempty',s.subs)), b=[]; return; end
      if ischar(s.subs{1}) && ~strcmp(s.subs{1},':'), b=get(b, s.subs{:}); return; end
      if length(s.subs) == 1 && all(s.subs{:} == 1), return; end
      d=get(b,'Signal'); d=d(s.subs{:});  b=set(b,'Signal', d);

      d=get(b,'Error');  if numel(d) > 1 & numel(d) == numel(get(a,'Error')), d=d(s.subs{:}); b=set(b,'Error', d); end

      d=get(b,'Monitor'); if numel(d) > 1 & numel(d) == numel(get(a,'Monitor')), d=d(s.subs{:}); b=set(b,'Monitor', d); end

      % must also affect axis
      for index=1:ndims(b)
        if index <= length(b.Alias.Axis)
          x = getaxis(b,index);
          ax= b.Alias.Axis{index};   % definition of Axis
          nd = size(x); nd=nd(find(nd>1));
          if length(size(x)) == length(size(b)) && ...
                 all(size(x) == size(b))  && all(length(nd) == length(s.subs)) % meshgrid type axes
            b = setaxis(b, index, ax, x(s.subs{:}));
          else  % vector type axes
            b = setaxis(b, index, ax, x(s.subs{index}));
          end
        end
      end 
      
      b = copyobj(b);
      
      % add command to history
      if ~isempty(inputname(2))
        toadd = [ inputname(2) ];
      elseif length(s.subs) == 1
        toadd = [  mat2str(double(s.subs{1})) ];
      elseif length(s.subs) == 2
        toadd = [  mat2str(double(s.subs{1})) ', ' mat2str(double(s.subs{2})) ];
      else
        toadd = [ '<not listable>' ];  
      end
      if ~isempty(inputname(1))
        toadd = [  b.Tag ' = ' inputname(1) '(' toadd ');' ];
      else
        toadd = [ b.Tag ' = ' a.Tag '(' toadd ');' ];
      end
  
      b = iData_private_history(b, toadd);
      % final check
      b = iData(b);
      % reset warnings
      try
        warning(warn.seta);
        warning(warn.geta);
        warning(warn.get);
      catch
        warning(warn);
      end

    end               % if single iData
  case '{}'
    if length(b(:)) > 1   % iData array
      b = b(s.subs{:});
    else
      if isnumeric(s.subs{:}) & length(s.subs{:}) == 1
        b=getaxis(b, s.subs{:});
      elseif ischar(s.subs{:}) & length(str2num(s.subs{:})) == 1
        b=getaxis(b, s.subs{:}); % definition of axis
      elseif ischar(s.subs{:})
        b=getalias(b, s.subs{:}); % same as b.'alias'
      else
        iData_private_error(mfilename, [ 'do not know how to extract cell index in ' inputname(1)  ' ' b.Tag '.' ]);
      end
    end
  case '.'
    if ~isstruct(b)
      b = get(b,s.subs);          % get field from iData
    else
      b = getfield(b,s.subs);     % get field from struct
    end
  end   % switch s.type
end % for s index level
