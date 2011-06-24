function b = subsasgn(a,S,val)
% b = subsasgn(a,index,b) : iData indexed assignement
%
%   @iData/subsasgn: function defines indexed assignement 
%   such as a(1:2,3) = b
%   The special syntax a{0} multiplies the value by the Monitor and then assigns 
%   the Signal, and a{n} assigns the axis of rank n.
%     When the assigned value is a char, the axis definition is set (as in setaxis).
%     When the assigned value is numeric, the axis value is set (as in set).
%   The special syntax a{'alias'} is a quick way to define an alias.
%
% Version: $Revision: 1.15 $
% See also iData, iData/subsref

% This implementation is very general, except for a few lines
% EF 27/07/00 creation
% EF 23/09/07 iData implementation

b = a;  % will be refined during the index level loop

if isempty(S)
  return
end

% first handle object array for first index
if length(b(:)) > 1 & (strcmp(S(1).type,'()') | strcmp(S(1).type,'{}'))
  c = b(S(1).subs{:}); % get elements in the array
  if ~isvector(c) d = c(:); else d=c; end
  if length(d) == 1
    d = iData(val);
  else
    if ~isempty( S(2:end) ) % array([1 2 3]).something = something
		  for j = 1:length(d)
		    d(j) = subsasgn(d(j),S(2:end),val);
		  end
		else	% single level array assigment as array([1 2 3]) = something
			for j = 1:length(d)
				if length(val) == 1, d(j) = iData(val);
				elseif length(val) == length(d)
					d = iData(val);
				else
					iData_private_error(mfilename, [ 'can not assign ' num2str(length(d)) ' iData array to ' num2str(length(val)) ' ' class(val) ' array for object ' inputname(1) ' ' b(1).Tag ]);
				end
			end
		end
  end
  if prod(size(c)) ~= 0 & prod(size(c)) == prod(size(d)) & length(c) > 1
    c = reshape(d, size(c));
  else
    c = d;
  end
  b(S(1).subs{:}) = c;
else
  i = 0;
  while i < length(S)     % can handle multiple index levels
    i = i+1;
    s = S(i);
    switch s.type
    case '()'       
      if length(b(:)) > 1   % array() -> deal on all elements
      % SYNTAX: array(index) = val: set Data using indexes
        c = b(:);           
        for j = 1:length(s.subs{:})
          c(j) = subsasgn(c(j),s,val);
        end
        b = reshape(c, size(b));
      elseif ~isa(val, 'iData') % single object() = Signal (val must be num)
      % SYNTAX: object(index) = numeric: set Signal, index can be a subset
        if ~isnumeric(val)
          iData_private_error(mfilename, [ 'object(' num2str(s.subs{:}) ') = ' class(val) ' but expects a numerical value to assign Signal in object ' inputname(1) ' ' b.Tag ]);
        end
        % this is where specific class structure is taken into account
        cmd=b.Command;
        d = get(b, 'Signal');
        d(s.subs{:}) = val;
        b = set(b, 'Signal', d);
        if isempty(val) % remove columns/rows in data: Update Error, Monitor and Axes
          for index=1:ndims(b)
            x = getaxis(b,index);
            if all(size(x) == size(b)) % meshgrid type axes
              x(s.subs{:}) = [];
              b = setaxis(b, index, x);
            else  % vector type axes
              x(s.subs{index}) = [];
            end
          end
        end
        
        % add command to history
        toadd = '(';
        if ~isempty(inputname(2))
          toadd = [ toadd inputname(2) ];
        elseif length(s.subs) == 1
          toadd = [ toadd mat2str(s.subs{1}) ];
        else
          toadd = [ toadd mat2str(s.subs{1}) ', ' mat2str(s.subs{2}) ];
        end
        toadd = [ toadd ') = ' ];
        if ~isempty(inputname(3))
          toadd = [ toadd inputname(3) ];
        elseif ischar(val)
          toadd = [ toadd '''' val '''' ];
        elseif isnumeric(val) | islogical(val)
          if length(size(val)) > 2, val=val(:); end
          if numel(val) > 10, 
            val=val(1:10); toadd = [ toadd mat2str(val) '...' ]; 
          else
            toadd = [ toadd mat2str(val) ];
          end
        else
          toadd = [ toadd  '<not listable>' ];
        end
        if ~isempty(inputname(1))
          toadd = [ inputname(1) toadd ';' ];
        else
          toadd = [ a.Tag toadd ';' ];
        end
        b.Command=cmd;
        b = iData_private_history(b, toadd);
        % final check
        b = iData(b);
      elseif length(s.subs{:}) == 1 && s.subs{:} == 1
      % SYNTAX: object(1) = iData: just assign objects
      	b = val;
      	return
      end                 % if single object
    case '{}'
      if length(b(:)) > 1   % object array -> deal on all elements
      % SYNTAX: array{ref}=val
        c = b(:);
        for j = 1:length(c)
          c(j) = subsasgn(c(j),s,val);
        end
        b = reshape(c, size(b));
      else
      % object{ref}
        if isnumeric(s.subs{:}) & length(s.subs{:}) == 1
        % object{axis_rank}
          if s.subs{:} == 0 & ischar(val)
          % SYNTAX: object{0} = numeric|char_link
            b = setalias(b, 'Signal', val);
            iData_private_warning(mfilename, [ 'Redefine Axis 0-th Signal in object ' inputname(1) ' ' b.Tag ]);
          else 
          % SYNTAX: object{axis_rank} = numeric|char
            if s.subs{:} <= length(b.Alias.Axis)
              if all(s.subs{:} > 0)
                ax= b.Alias.Axis{s.subs{:}}; % re-definition of Axis
              else
                ax='Signal';
              end
              if isempty(ax) & isnumeric(val) % change numerical value of axis
              % SYNTAX: object{axis_rank} = numeric
                ax=[ 'Axis_' num2str(num2str(s.subs{:})) ];
                % need to create this axis
                setalias(b, ax, val);
                setaxis(b, s.subs{:}, ax);
              else  % change axis definition
              % SYNTAX: object{axis_rank} = char: redefine axis link
              % SYNTAX: object{axis_rank} = other: redefine alias
                if ischar(val), b = setaxis(b, s.subs{:}, val);
                else 
                  % special case for Signal, which should take into account the Monitor
                  if strcmp(ax, 'Signal')
                    m  = get(b,'Monitor'); m=real(m);
                    if not(all(m == 1 | m == 0))
                      val = genop(@times,val,m);
                    end
                  end
                  b = set(b, ax, val); 
                end
              end
            else 
            % create new axis as axis rank exceeds dim(object)
              if isnumeric(val)
              % SYNTAX: object{axis_rank > dim} = numeric
                iData_private_warning(mfilename, [ num2str(s.subs{:}) '-th rank Axis  has not been defined yet. Defining "Axis_' num2str(num2str(s.subs{:})) '" in object ' inputname(1) ' ' b.Tag ]);
                ax=[ 'Axis_' num2str(num2str(s.subs{:})) ];
                setalias(b, ax, val);
                setaxis(b, s.subs{:}, ax);
              elseif ischar(val)
              % SYNTAX: object{axis_rank > dim} = char
                ax=s.subs{:}; 
              else
                iData_private_error(mfilename, [ num2str(s.subs{:}) '-th rank Axis can not be assigned in object ' inputname(1) ' ' b.Tag ]);
              end
              if ischar(val), b = setaxis(b, s.subs{:}, val);
              else b = set(b, ax, val); end
            end
          end
        elseif ischar(s.subs{:}) & length(str2num(s.subs{:}))
        % SYNTAX: object{'axis_rank'} = val -> object{axis_rank} = val
          b=setaxis(b, str2num(s.subs{:}), val);
        elseif ischar(s.subs{:})
        % SYNTAX: object{'field'} = val
          b=setalias(b, s.subs{:}, val);
        else
        % SYNTAX: object{index} = other: will probably fail...
          b{s.subs{:}} = val;       
        end
      end
    case '.'
    % SYNTAX: object.field = val
      if length(b(:)) > 1   % object array -> deal on all elements
        c = b(:);
        for j = 1:length(c)
          c(j) = subsasgn(c(j),s,val);
        end
        b = reshape(c, size(b));
      else
      % ensure recursive assignment (see Matlab doc)
        if i < length(S)
          next_s = S(i+1);
          if strcmp(next_s.type, '()') | strcmp(next_s.type, '{}')
            tmp = get(b,s.subs);
            tmp(next_s.subs{:}) = val;
            b = set(b, s.subs, tmp);
            i = i + 2;  % jump next
          else
            if strcmp(next_s.type, '.')
              c = getfield(b, s.subs);
              c = setfield(c, next_s.subs, val);
              b = setfield(b, s.subs, c);
              i = i+1;
            else
              iData_private_error(mfilename, [ 'can not handle ' next_s.type ' type subscript for object ' inputname(1) ' ' b.Tag ]);
            end
          end
        else
          b = set(b,s.subs,val);    % set field
        end
      end
    end   % switch s.type
  end % while s index level
  
end

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),b);
end
