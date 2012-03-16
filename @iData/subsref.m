function b = subsref(a,S)
% b = subsref(a,s) : iData indexed references
%
%   @iData/subsref: function returns subset indexed references
%     such as a(1:2) or a.field.
%   The special syntax a{0} where a is a single iData returns the 
%     Signal/Monitor, and a{n} returns the axis of rank n.
%
% Version: $Revision: 1.28 $
% See also iData, iData/subsasgn

% This implementation is very general, except for a few lines
% EF 27/07/00 creation
% EF 23/09/07 iData impementation
% ==============================================================================
% inline: private function iData_getAliasValue (mainly used)
% calls:  subsref (recursive), getaxis, getalias, get(Signal, Error, Monitor)

b = a;  % will be refined during the index level loop

if isempty(S)
  return
end

for i = 1:length(S)     % can handle multiple index levels
  s = S(i);
  switch s.type
  case '()' % ======================================================== array
    if numel(b) > 1   % iData array
      b = b(s.subs{:});
    else                  % single iData
      % this is where specific class structure is taken into account
      if any(cellfun('isempty',s.subs)), b=[]; return; end        % b([])
      if ischar(s.subs{1}) && ~strcmp(s.subs{1},':')              % b(name) -> s.(name) alias/field value
        s.type='.';
        b=subsref(b, s); return;
      end
      if length(s.subs) == 1 && all(s.subs{:} == 1), continue; end  % b(1)
      
      iData_private_warning('enter',mfilename);
      
      ds=iData_getAliasValue(b,'Signal'); 
      d=ds(s.subs{:});                          % b(indices)
      b=set(b,'Signal', d);  b=setalias(b,'Signal', d);
      clear ds
      
      de=iData_getAliasValue(b,'Error'); 
      if numel(de) > 1 && isnumeric(de) 
        try % in case Error=sqrt(Signal), the Error is automatically changed when Signal is -> fail
          d=de(s.subs{:}); b=set(b,'Error', d); b = setalias(b, 'Error', d);
        end
      end
      clear de

      dm=iData_getAliasValue(b,'Monitor');
      if numel(dm) > 1 && isnumeric(dm)
        d=dm(s.subs{:}); b=set(b,'Monitor', d);  b = setalias(b, 'Monitor', d);
      end
      clear dm

      % must also affect axis
      for index=1:ndims(b)
        if index <= length(b.Alias.Axis)
          x = getaxis(b,index);
          ax= b.Alias.Axis{index};   % definition of Axis
          nd = size(x); nd=nd(nd>1);
          if length(size(x)) == length(size(b)) && ...
                 all(size(x) == size(b))  && all(length(nd) == length(s.subs)) % meshgrid type axes
            b = setaxis(b, index, ax, x(s.subs{:}));
          elseif max(s.subs{index}) <= numel(x) % vector type axes
            b = setaxis(b, index, ax, x(s.subs{index}));
          else
            iData_private_warning(mfilename,[ 'The Axis ' num2str(size(index)) ' [' ...
    num2str(size(x)) ' can not be resize as a [' num2str(size(s.subs{index})) ...
    '] vector in iData object ' b.Tag ' "' b.Title '".\n\tTo use the default Error=sqrt(Signal) assign s.Error=[].' ]);
          end
        end
      end 
      
      b = copyobj(b);
      
      % add command to history
      if ~isempty(inputname(2))
        toadd =   inputname(2);
      elseif length(s.subs) == 1
        toadd =    mat2str(double(s.subs{1}));
      elseif length(s.subs) == 2
        toadd = [  mat2str(double(s.subs{1})) ', ' mat2str(double(s.subs{2})) ];
      else
        toadd =   '<not listable>';  
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
      iData_private_warning('exit',mfilename);

    end               % if single iData
  case '{}' % ======================================================== cell
    if isnumeric(s.subs{:}) && isscalar(s.subs{:})
      b=getaxis(b, s.subs{:});  % b{rank} value of axis
    elseif ischar(s.subs{:}) && ~isnan(str2double(s.subs{:}))
      b=getaxis(b, s.subs{:});  % b{'rank'} definition of axis
    elseif ischar(s.sub{:})
      b=getalias(b, s.subs{:}); % b{'alias'} same as b.'alias' definition
    else
      iData_private_error(mfilename, [ 'do not know how to extract cell index in ' inputname(1)  ' ' b.Tag '.' ]);
    end
  case '.'  % ======================================================== structure
    % protect some fields
    fieldname = s.subs;
    if length(fieldname) > 1 && iscell(fieldname)
      fieldname = fieldname{1};
    end
    if strcmpi(fieldname, 'filename') % 'alias of alias'
      fieldname = 'Source';
    elseif strcmpi(fieldname, 'history')
      fieldname = 'Command';
    elseif strcmpi(fieldname, 'axes')
      fieldname = 'Axis';
    end
    if any(strcmpi(fieldname, 'alias'))
      b = getalias(b);
    elseif any(strcmpi(fieldname, 'axis'))
      b = getaxis(b);
    elseif any(strcmpi(fieldname, fieldnames(b))) % structure/class def fields: b.field
      b = b.(fieldname);
      if isnumeric(b) && any(strcmp(fieldname, {'Date','ModificationDate'}))
        b = datestr(b);
      end
    else
      b = iData_getAliasValue(b,fieldname); % get alias value from iData: b.alias MAIN SPENT TIME
    end
    
    % test if the result is again an Alias or Field
    if ischar(b) && size(b,1) == 1
      if any(strcmpi(b, fieldnames(a(1))))
        b = a.(b);      % fast access to static fields
      elseif any(strcmpi(strtok(b,'.'), fieldnames(a))) || any(strcmpi(strtok(b,'.'), a(1).Alias.Names))
        b = get(a, b);  % try to evaluate char result
      end
    end
  end   % switch s.type
end % for s index level

% ==============================================================================
% private function iData_getAliasValue
function val = iData_getAliasValue(this,fieldname)
% iData_getAliasValue: iData alias evaluation (not the link, but the value)
%   evaluates s.name to be first s.link, then 'link' (with 'this' defined).
%   NOTE: for standard Aliases (Error, Monitor), makes a dimension check on Signal

% EF 23/09/07 iData impementation
  val = [];
  if ~isa(this, 'iData'),   return; end
  if ~isvarname(fieldname), return; end % not a single identifier (should never happen)

  % searches if this is an alias (it should be)
  alias_num   = find(strcmpi(fieldname, this.Alias.Names));  % index of the Alias requested
  if isempty(alias_num), 
    iData_private_error(mfilename, sprintf('can not find Property "%s" in object %s "%s".', fieldname, this.Tag, this.Title ));
    return; 
  end                    % not a valid alias
  
  alias_num = alias_num(1);
  name      = this.Alias.Names{alias_num};
  val       = this.Alias.Values{alias_num};  % definition/value of the Alias

  if (~isnumeric(val) && ~islogical(val))
    % the link evaluation must be numeric in the end...
    if ~ischar(val),       return; end  % returns numeric/struct/cell ... content as is.
    if  strcmp(val, name), return; end  % avoids endless iteration.
    
    % val is now only a char
    % handle URL content (possibly with # anchor)
    if  (strncmp(val, 'http://', length('http://'))  || ...
         strncmp(val, 'https://',length('https://')) || ...
         strncmp(val, 'ftp://',  length('ftp://'))   || ...
         strncmp(val, 'file://', length('file://')) )
      % evaluate external link
      val = iLoad(val); % stored as a structure
      return
    end
    
    % gets the alias value (evaluate the definition) this.alias -> this.val
    if ~isempty(val)
      % handle # anchor style alias
      if val(1) == '#', val = val(2:end); end % HTML style link
      % evaluate the alias definition (recursive call through get -> subsref)
      try
        % in case this is an other alias/link: this is were we propagate in the object
        val = get(this,val); % gets this.(val)                    MAIN SPENT TIME
      catch
        % evaluation failed, the value is the char (above 'get' will then issue a

        % 'can not find Property' error, which will come there in the end
      end
    end
  end

  % link value has been evaluated, do check in case of standard aliases
  if strcmp(fieldname, 'Error')         % Error is sqrt(Signal) if not defined 
    if ~isempty(val)
      if all(val(:) == val(end))
        val = val(end);
      end
    else
      s = iData_getAliasValue(this,'Signal');
      if isnumeric(s)
        val = sqrt(abs(double(s)));
      end
    end
    if ~isempty(val) && ~isscalar(val) && ~isequal(size(val),size(this))
      iData_private_warning(mfilename,[ 'The Error [' num2str(size(val)) ...
      '] has not the same size as the Signal [' num2str(size(this)) ...
      '] in iData object ' this.Tag ' "' this.Title '".\n\tTo use the default Error=sqrt(Signal) assign s.Error=[].' ]);
    end
  elseif strcmp(fieldname, 'Monitor')  % Monitor is 1 by default
    if isempty(val), val=1;
    elseif all(val(:) == val(end))
      val = val(end);
    end
    if val == 0, val=1; end
    if ~isempty(val) && length(val) ~= 1 && ~all(size(val) == size(this))
      iData_private_warning(mfilename,[ 'The Monitor [' num2str(size(val)) ...
        '] has not the same size as the Signal [' num2str(size(this)) ...
        '] in iData object ' this.Tag ' "' this.Title '".\n\tTo use the default Monitor=1 use s.Monitor=[].' ]);
    end
  end

