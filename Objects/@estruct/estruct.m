classdef estruct < dynamicprops
% ESTRUCT Create or convert to an estruct object.
%   S = ESTRUCT('field1',VALUES1,'field2',VALUES2,...) creates
%   an object with the specified fields and values (as properties).
%   ESTRUCT is similar to STRUCT, but is designed to hold scientific data.
%   After creation, it can be manipulated as a normal array.
% 
%   ESTRUCT(a,b,c,...) imports input arguments into an estruct object.
%   ESTRUCT(X,Y,... Signal) sets the Signal to the last given numeric array, and
%   the axes as the previous vectors/arrays.
% 
%   ESTRUCT('filename') imports the file name and store its content into a Data property.
%   Multiple files can be given, each producing a distinct object arranged
%   into an array. File names can also be given as URL (file: http: https: ftp:)
%   and/or be compressed (zip, gz, tar, Z). The import uses a guessed importer
%   and does not apply any post-process filter (raw data from file).
%   Use 'LOAD(estruct, file, loader)' to specify the importer and apply
%   post-processing.
%
%   ESTRUCT silent|verbose|debug sets the verbosity level. The verbosity level is
%   common to all objects and methods. You may also change the verbosity (level
%   of output messages) with the 'verbose' property, e.g.
%     a=estruct; a.verbose = 1;
%   A null verbosity corresponds with no messages (silent). 
%
%   ESTRUCT(OBJ) when OBJ is an ESTRUCT checks for object consistency and return.
%   ESTRUCT(OBJ) converts the object OBJ into its equivalent
%   estruct.  The initial class information is lost.
% 
%   ESTRUCT([]) creates an empty object.
%
%   Example: s = estruct('type',{'big','little'},'color','red','x',{3 4}); isstruct(s)
%   Example: s = estruct(1:10); isstruct(s)
%   Version: $Date$ $Version$ $Author$
%   See also estruct/load, isstruct, setfield, getfield, fieldnames, orderfields,
%     isfield, rmfield, deal, substruct, struct2cell, cell2struct.

properties

    % MetaData properties (sorted in alpha order)
    Creator         = mfilename;  % The Creator/Software of the object
    Command         ={};          % The history of the object (cellstr)
    Date            = clock;      % The Date of creation of the object
    Data            =[];          % Where we store most of the Data
    Label           ='';          % The label of the object
    ModificationDate=[];          % The Date of last modification
    Name            ='';          % The Name/Title of the object
    Source          = pwd;        % The data source
    User            = getenv('USER'); % Who created the object
    UserData        =[];          % An area to store what you want
    Tag             ='';          % A unique ID for the object

  end % properties

  properties (Access=private, Hidden=true)     % internal use

    Private
    % Data handling: Signal, Axes, ...
    Labels  = struct(); % struct: Labels.Signal, ... Labels.Axes{1:ndims}
    Axes    = {};       % cell{1:ndims} e.g. aliases
    verbose = 0;        % 0:silent, 1:normal, 2:debug
    Attributes = [];    % field attributes. Map the object structure.
    DisplayName='';     % The apparent name of the object, e.g. for plots. Alias to Label.
  end

  properties (Access=protected, Constant=true)  % shared by all instances
    properties_Protected={'properties_Protected','properties_Base', ...
      'Axes','Tag','Private','Labels' } % can not be changed manually
    properties_Base={'Creator', 'Command', 'Date', 'Data', 'DisplayName', ...
      'Label', 'ModificationDate', 'Source', 'Tag', 'Name', 'User', 'UserData', ...
      'Attributes' };
  end
  
  events
    % These are inherited from dynamicprops
    % PropertyAdded
    % PropertyRemoved
    % ObjectBeingDestroyed
    ObjectUpdated % triggered when a check is done
    ObjectPlotted % triggered when the object is plotted
  end % events

% ------------------------------------------------------------------------------

  methods
    function new = estruct(varargin)
    % ESTRUCT Create or convert to an estruct object.
    %   S = ESTRUCT('field1',VALUES1,'field2',VALUES2,...) creates
    %   an object with the specified fields and values (as properties).
    %   ESTRUCT is similar to STRUCT, but is designed to hold scientific data.
    %   After creation, it can be manipulated as a normal array.
    % 
    %   ESTRUCT(a,b,c,...) imports input arguments into an estruct object.
    %   ESTRUCT(X,Y,... Signal) sets the Signal to the last given numeric array, and
    %   the axes as the previous vectors/arrays.
    % 
    %   ESTRUCT('filename') imports the file name and store its content into a Data property.
    %   Multiple files can be given, each producing a distinct object arranged
    %   into an array. File names can also be given as URL (file: http: https: ftp:)
    %   and/or be compressed (zip, gz, tar, Z). The import uses a guessed importer
    %   and does not apply any post-process filter (raw data from file).
    %   Use 'LOAD(estruct, file, loader)' to specify the importer and apply
    %   post-processing.
    %
    %   ESTRUCT silent|verbose|debug sets the verbosity level. The verbosity level is
    %   common to all objects and methods. You may also change the verbosity (level
    %   of output messages) with the 'verbose' property, e.g.
    %     a=estruct; a.verbose = 1;
    %   A null verbosity corresponds with no messages (silent). 
    %
    %   ESTRUCT(OBJ) when OBJ is an ESTRUCT checks for object consistency and return.
    %   ESTRUCT(OBJ) converts the object OBJ into its equivalent
    %   estruct.  The initial class information is lost.
    % 
    %   ESTRUCT([]) creates an empty object.
    %
    % Example: s = estruct('type',{'big','little'},'color','red','x',{3 4}); ...
    %          isstruct(s)
    % Version: $Date$ $Version$ $Author$
    % See also estruct.load, isstruct, setfield, getfield, fieldnames, orderfields,
    %   isfield, rmfield, deal, substruct, struct2cell, cell2struct.

      persistent meth

      if isempty(meth), meth = methods(mfilename); end
      
      if nargin == 1 && isa(varargin{1}, 'estruct')
        new = axescheck(varargin{1}, 'force');
        return
      end

      warning('off','MATLAB:structOnObject');
      new.Private.cache = []; % init cache to empty
      new.Private.cache.methods = meth;
      % handle Tag number
      new.Tag = [ 'iD' sprintf('%0.f', private_id()) ]; % unique ID
      new.ModificationDate = new.Date;

      % add our 'static' properties so that they are equally handled by
      % subsref/subsasgn
      new.addprop('Error');           % e.g. an alias or empty
      new.addprop('Monitor');         % e.g. an alias or empty
      new.addprop('Signal');          % e.g. an alias
      new.Error = 'matlab: sqrt(this.Signal)';
      if ~nargin, new.Private.cache.check_requested = false; return; end

      % collect items to store: as structures, as data files, and others
      structs = {}; % cell: will contain struct('name','value')
      flag_num_arg = 0;

      % append arguments as properties (not from files)
      index=1;
      while index<=nargin
        this = varargin{index}; s = [];
        if isobject(this)
          if numel(this) == 1, this = struct(this); else this = []; end
        end
        if ischar(this) && ~isempty(this)
          switch lower(this)
          case {'verbose','normal','default'}
            new.verbose = 1;
            if exist('iLoad','file'), iLoad('verbose'); end
            new.Private.cache.check_requested = false;
            return;
          case 'debug'
            new.verbose = 2;
            if exist('iLoad','file'), iLoad('debug'); end
            new.Private.cache.check_requested = false;
            return;
          case 'silent'
            new.verbose = 0;
            if exist('iLoad','file'), iLoad('silent'); end
            new.Private.cache.check_requested = false;
            return;
          end

          if index<numel(varargin) && isvarname(this)   % input: name/value pair, e.g. 'par',value, ...
            s.name = this;
            s.value= varargin{index+1};
            structs{end+1} = s; this = [];
            varargin{index} = []; varargin{index+1} = []; % clear memory
            index=index+1;  % increment name/value pair
          elseif ~isempty(dir(this)) || any(strcmp(strtok(this, ':'), {'http' 'https' 'ftp' 'file'}))
            % pass: we handle this below when assembling the object(s)
            this = [];
          else
            s.value = str2struct(this); % convert to struct ?
            s.name  = inputname(index);
            structs{end+1} = s; this = [];
            varargin{index} = [];
          end
        elseif isstruct(this) && isfield(this, 'Source') && isfield(this, 'Loader') % iLoad struct
          % pass: we handle this iLoad struct when assembling the object(s)
            this = [];
        end

        if ~isempty(this)
          s.value = this;
          s.name  = inputname(index);
          if isnumeric(this) | islogical(this)
            flag_num_arg = flag_num_arg+1;
            s.axis=flag_num_arg;
          end
          structs{end+1} = s;
        end

        index=index+1;
      end % varargin

      % fill in direct name/value pairs stored into 'structs'
      for index=1:numel(structs)
        s = structs{index};
        if isempty(s.name), s.name=sprintf('%s_%i', class(s.value), index); end
        if ~isfield(new, s.name), new.addprop(s.name); end
        new.(s.name)=s.value;
        if (isnumeric(s.value) | islogical(s.value)) && ~isscalar(s.value)
          history(new, 'set', new, s.name, s.value);
          if isfield(s,'axis')
            if s.axis == flag_num_arg && ~strcmp(s.name,'Signal') % last num -> Signal
              s.Signal = s.name; % alias
            else
              s.Axes{s.axis} = s.name;
            end
          end
        end
        structs{index} = []; % clear memory
      end
      numel_new = 1;
      new0 = copyobj(new);

      % now build the final 'master' object
      for index_varg=1:numel(varargin) % loop on initial input arguments for iLoad
        arg = varargin{index_varg};
        if isempty(arg) || (~ischar(arg) && ~iscellstr(arg) && ~isstruct(arg)), continue; end
        if ischar(arg), arg = cellstr(arg); end
        for index_arg=1:numel(arg)  % loop on input argument content when e.g. cellstr array
          if iscell(arg)
            this = arg{index_arg}; % must be a char or single cellstr
          else
            this = arg(index_arg); % e.g. struct/object
          end
          if isempty(this), continue; end
          if ischar(this) && (~isempty(dir(this)) || any(strcmp(strtok(this, ':'), {'http' 'https' 'ftp' 'file'})))
            try
              this = load(estruct, this, 'raw'); % import file/URL, no post-processing
            catch ME
              disp([ mfilename ': WARNING: failed importing ' char(this) ]);
            end
          end
          % 'this' can be initial char/cellstr, a struct from iLoad, or a cell of structs from iLoad
          if ~iscell(this), this = { this }; end
          for index_data=1:numel(this)
            if numel(this{index_data}) == 0, continue; end
            if isa(this{index_data}, 'estruct')
              new1 = this{index_data};
            elseif isempty(this{index_data}), continue;
            else
              if numel_new == 1
                new1 = new;
              else
                new1 = copyobj(new0);
              end

              if isstruct(this{index_data})
                new1 = struct2estruct(this{index_data}, new1);
              else
                new1.Data = this{index_data}; 
              end
              history(new1, mfilename, this{index_data});
              this{index_data} = []; % and clear memory
              if isstruct(new1.Data) && isfield(new1.Data, 'Source') && isfield(new1.Data, 'Loader') % iLoad struct
                % we transfer the iLoad struct to the estruct.
                new1 = struct2estruct(new1.Data, new1); % updates estruct (in 'private')
              end
            end % isa estruct or other type to concatenate
            
            if numel_new == 1
              new = new1;
            else
              new = [ new ; new1(:) ]; % build array
            end % index_data
            numel_new = numel_new+1;
          end
        end % index_arg (content of varg in case this is an array/cellstr for iLoad)

      end % index_varg (initial input arg)

      % final check
      for index_new = 1:numel(new) % post process may create more objects
        set(new(index_new), 'Private.cache.check_requested',true); % request a check at first 'get'
      end

    end % estruct (instantiate)

% ------------------------------------------------------------------------------

    function tf = isstruct(self)
    %  ISSTRUCT True for structures.
    %      ISSTRUCT(S) returns logical true (1) if S is a structure
    %      and logical false (0) otherwise.
    %
    % Example: s = estruct; isstruct(s)
    %
    % See also estruct, isfield, iscell, isnumeric, isobject.
      tf = true; % isa(self, mfilename); always true
    end

    function self=rmfield(self, f)
    %   RMFIELD Remove fields from a structure array.
    %      S = RMFIELD(S,'field') removes the specified field from the
    %      m x n structure array S. The size of input S is preserved.
    %
    %      S = RMFIELD(S,FIELDS) removes more than one field at a time
    %      when FIELDS is a character array or cell array of strings.  The
    %      changed structure is returned. The size of input S is preserved.
    %
    %      S = RMFIELD(S,'all') removes all fields except base ones. The resulting
    %      object retains metadata, but removes any additional alias/property.
    %
    % Example: s=estruct; s.test =1; ...
    %          isfield(s,'test') && ~isfield(rmfield(s,'test'),'test')
    %
    %      See also setfield, getfield, isfield, fieldnames.
      if nargin ~= 2, return; end
      if ischar(f) && strcmp(f, 'all')
        f = fieldnames(self);
      end
      if numel(self) == 1
        if ~iscellstr(f), f = cellstr(f); end
        for index=1:numel(f)
          if isfield(self, f{index})
            if ~any(strcmp(f{index}, self.properties_Base)) ...
            && ~any(strcmp(f{index}, self.properties_Protected)) ...
            && ~any(strcmp(f{index}, {'Signal','Error','Monitor'}))
              delete(findprop(self, f{index}));
            end
          end
        end
      else
        for index=1:numel(self); rmfield(self(index), f); end
      end
    end % rmfield

    function self=rmalias(self, varargin)
    %   RMALIAS Remove fields from a structure array.
    %      S = RMALIAS(S,'field') removes the specified field from the
    %      m x n structure array S. The size of input S is preserved.
    %
    %      S = RMALIAS(S,FIELDS) removes more than one field at a time
    %      when FIELDS is a character array or cell array of strings.  The
    %      changed structure is returned. The size of input S is preserved.
    %
    %      S = RMALIAS(S,'all') removes all fields except base ones. The resulting
    %      object retains metadata, but removes any additional alias/property.
    %
    % Example: s=estruct; s.test =1; ...
    %          isfield(s,'test') && ~isfield(rmalias(s,'test'),'test')
    %
    %      See also setalias, getalias, isfield, fieldnames.

    % compatibility with original iData (2007-2019)
      self = rmfield(self, varargin{:});
    end % rmalias

    function self=rmaxis(self, f)
    %   RMAXIS Remove an axis from object(s).
    %     RMAXIS(S, AX) removes the axis AX specified as a single rank index.
    %
    % Example: s=estruct(1:10); s{1}=2:11; rmaxis(s,1); all(s{1} == 1:10)
      if nargin ~= 2 || ~isnumeric(f), return; end
      if numel(self) == 1
        f = find( 1 <= f & f <= numel(self.Axes) );
        for index=1:numel(f); self.Axes{f(index)} = []; end
      else
        for index=1:numel(self); rmaxis(self(index), f); end
      end
    end % rmaxis

    function c = struct2cell(self)
    %   STRUCT2CELL Convert structure array to cell array.
    %      C = STRUCT2CELL(S) converts the M-by-N structure S (with P fields)
    %      into a P-by-M-by-N cell array C.
    %
    %      If S is N-D, C will have size [P SIZE(S)].
    %
    %      Example: c = struct2cell(estruct); ...
    %               iscell(c)
    %        clear s, s.category = 'tree'; s.height = 37.4; s.name = 'birch';
    %        c = struct2cell(s); f = fieldnames(s);
    %
    %      See also cell2struct, fieldnames.

      if numel(self) == 1
        c = struct2cell(struct(self));
      else
        c = arrayfun('struct2cell', self);
      end

    end % struct2cell

    function new = cell2struct(self, varargin)
    %   CELL2STRUCT Convert cell array to estruct array.
    %      S = CELL2STRUCT(C,FIELDS,DIM) converts the cell array C into
    %      the structure S by folding the dimension DIM of C into fields of
    %      S.  SIZE(C,DIM) must match the number of field names in FIELDS.
    %      FIELDS can be a character array or a cell array of strings.
    %
    %      Example: s = cell2struct(estruct, {'tree',37.4,'birch'}, ...
    %         {'category','height','name'},2); isa(s, 'estruct')
    %
    %      See also struct2cell, fieldnames.
      if nargin < 3, d=2; end
      if nargin < 2, error([ mfilename ': cell2struct: requires at least 2 arguments' ]); end
      if numel(self) > 1
        error([ mfilename ': cell2struct(estruct, ...) only works with a single estruct.' ])
      end
      new = copyobj(self, cell2struct(varargin{:}));
    end % cell2struct

    function s = cellstr(self, varargin)
      %  CELLSTR Convert an estruct object into a cell array of strings.
      %      S = CELLSTR(X) converts the structure X into a cell of strings.
      %
      %      S = CELLSTR(X,'short') converts the structure X into a cellstr
      %      and a compact storage. The resulting representation does not
      %      allow to rebuild the full object.
      %
      %  Example: iscellstr(cellstr(estruct(1:10)))
      %
      % See also: char
      s = copyobj(self); s.Private = [];
      s = class2str(s, varargin{:});
      s = textscan(s, '%s','Delimiter',';');
      if numel(s) == 1, s=s{1}; end
    end % cellstr

    function v = cast(self, typ)
      %   CAST  Cast object to a different data type or class.
      %     B = CAST(A,NEWCLASS) casts A to class NEWCLASS. A must be convertible to
      %     class NEWCLASS. NEWCLASS must be the name of one of the builtin data types.
      %
      % Example: s=estruct(1:10); isa(cast(s, 'uint64'), 'uint64')
      %
      % See also: estruct.double, estruct.int32, estruct.uint32
      if nargin < 2, typ=''; end
      if isempty(typ), typ='double'; end
      if numel(self) == 1
        v = cast(subsref(self,struct('type','.','subs','Signal')), typ);
      else
        v = cell(size(self));
        for index=1:numel(self); v{index} = cast(self(index), typ); end
      end
    end % cast

    function v = double(self)
      %   DOUBLE Convert to double precision.
      %      DOUBLE(X) returns the double precision value of the object.
      %
      % Example: s=estruct(logical(-5:5)); isa(double(s), 'double')
      %
      % See also: estruct.single, estruct.cast
      v = cast(self, 'double');
    end

    function v = single(self)
      %   SINGLE Convert to SINGLE precision.
      %      SINGLE(X) returns the single precision value of the object.
      %
      % Example: s=estruct(logical(-5:5)); isa(single(s), 'single')
      %
      % See also: estruct.double, estruct.cast
      v = cast(self, 'single');
    end

    function v = logical(self)
      %   LOGICAL Convert to logical (0=false, any other=true).
      %      LOGICAL(X) returns the logical value of the object.
      %
      % Example: s=estruct(-5:5); isa(logical(s), 'logical')
      %
      % See also: estruct.double, estruct.cast
      v = cast(self, 'logical');
    end

    function v = uint32(self)
      %   UINT32 Convert to uint32, e.g. for indexing.
      %      UINT32(X) returns the unsigned integer alue of the object.
      %
      % Example: s=estruct(-5:5); isa(uint32(s), 'uint32')
      %
      % See also: estruct.double, estruct.cast
      v = cast(self, 'uint32');
    end

    function v = int32(self)
      %   INT32 Convert to int32.
      %      INT32(X) returns the integer alue of the object.
      %
      % Example: s=estruct(-5:5); isa(int32(s), 'int32')
      %
      % See also: estruct.double, estruct.cast
      v = cast(self, 'int32');
    end

    function s = repmat(self, varargin)
      % REPMAT Replicate and tile an array.
      %    B = repmat(A,M,N) creates a large matrix B consisting of an M-by-N
      %    tiling of copies of A. The size of B is [size(A,1)*M, size(A,2)*N].
      %    The statement repmat(A,N) creates an N-by-N tiling.
      %
      %    B = REPMAT(A,[M N]) accomplishes the same result as repmat(A,M,N).
      %
      %    B = REPMAT(A,[M N P ...]) tiles the array A to produce a
      %    multidimensional array B composed of copies of A. The size of B is
      %    [size(A,1)*M, size(A,2)*N, size(A,3)*P, ...].
      %
      % Example: numel(repmat(estruct, 11,10)) == 110
      %
      % See also: estruct.reshape, estruct.ones
      if numel(self) > 1
        error([ mfilename ': repmat(estruct, M,N,...) only works with a single estruct.' ])
      end
      % we create M*N copies and then reshape the array.
      s = [];
      for index=1:prod([ varargin{:}])
        this = copyobj(self);
        if index == 1, s = this; else s = [ s this]; end
      end
      s = reshape(s, varargin{:});
    end % repmat

    function o = ones(self, varargin)
      % ONES   Ones array.
      %    ONES(S,N) is an N-by-N matrix of 'S' structure.
      %
      %    ONES(S, M,N) or ONES(S, [M,N]) is an M-by-N matrix of 'S' structure.
      %
      %    ONES(S,M,N,P,...) or ONES(S,[M N P ...]) is an M-by-N-by-P-by-... array of
      %    'S' structure.
      %
      % Example: numel(ones(estruct, 11,10)) == 110
      %
      % See also: estruct.reshape
      if nargin == 1, o=self(1); return; end
      o = repmat(self, varargin{:});
    end % ones

    function z = zeros(self, varargin)
      % ZEROS  Zeros array.
      %    ZEROS(S,N) is an N-by-N matrix of empty structures.
      %
      %    ZEROS(S,M,N) or ZEROS(S,[M,N]) is an M-by-N matrix of empty structures.
      %
      %    ZEROS(S,M,N,P,...) or ZEROS(S,[M N P ...]) is an M-by-N-by-P-by-... array of
      %    empty structures.
      %
      %    ZEROS(S) removes all properties except base ones, and keep metadata.
      %
      % Example: numel(zeros(estruct, 11,10)) == 110
      %
      % See also: estruct.reshape
      if nargin == 1,
        z=estruct;
        z.Creator         = self.Creator;
        z.Date            = self.Date;
        z.DisplayName     = self.DisplayName;
        z.Source          = self.Source;
        z.Name            = self.Name;
        z.User            = self.User;
        z.UserData        = self.UserData;
        return
      end
      z = ones(estruct, varargin{:});
    end % zeros

    function s = struct(self)
      % STRUCT Create or convert to structure array.
      %   STRUCT(OBJ) converts the object OBJ into its equivalent
      %   structure.  The class information is lost.
      %
      % Example: s=estruct(1:10); isstruct(struct(s))
      if numel(self) == 1
        s = builtin('struct',self);
      else
        s = arrayfun('struct', self);
      end
    end % struct

    function tf = ismethod(self, m)
      % ISMETHOD  True if method of object.
      %   ISMETHOD(OBJ,NAME) returns 1 if string NAME is a method of object
      %   OBJ, and 0 otherwise.
      %
      % Example: s=estruct(1:10); ismethod(s, 'ismethod')
      try
        tf = any(strcmp(m ,self.Private.cache.methods));
      catch
        tf = any(strcmp(m, methods(self)));
      end
    end % ismethod

    % special SETTERS and GETTERS ----------------------------------------------
    
    function self = set.verbose(self, value)
    % SET VERBOSE Set the verbosity level
    %     0 silent
    %     1 normal
    %     2 verbose
    %     3 debug

      self.verbose = private_verbose(value); % shared verbose property
    end
    
    function value = get.verbose(self)
    % GET VERBOSE Get the verbosity level
    %     0 silent
    %     1 normal
    %     2 verbose
    %     3 debug

      value = private_verbose;                % shared verbose property
    end

    function self = set.DisplayName(self, value)
      self.Label = value; % alias to label property
    end

    function value = get.DisplayName(self)
      value = self.Label; % alias to label property
    end
    

  end % methods

% ------------------------------------------------------------------------------
end % classdef
