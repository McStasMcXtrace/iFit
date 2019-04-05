classdef estruct < dynamicprops
%ESTRUCT Create or convert to estruct.
%  S = ESTRUCT('field1',VALUES1,'field2',VALUES2,...) creates
%  an object with the specified fields and values (as properties).
%
%  ESTRUCT(a,b,c,...) imports input arguments into an estruct object.
%
%  ESTRUCT('filename') imports the file name and store its content into a Data property. 
%  Multiple files can be given, each producing a distinct object arranged
%  into an array. File names can also be given as URL (file: http: https: ftp:)
%  and/or be compressed (zip, gz, tar, Z).
%
%  ESTRUCT(OBJ) converts the object OBJ into its equivalent
%  estruct.  The initial class information is lost.
%
%  ESTRUCT([]) creates an empty object.
%
%  ESTRUCT is similar to STRUCT, but is designed to hold scientifi data.
%
% Example: s = estruct('type',{'big','little'},'color','red','x',{3 4}); isstruct(s)
% Version: $Date$ $
% See also isstruct, setfield, getfield, fieldnames, orderfields, 
%   isfield, rmfield, deal, substruct, struct2cell, cell2struct.

properties

    % MetaData properties (sorted in alpha order)
    Creator         = mfilename;
    Command         ={};
    Date            = clock;
    Data            =[]; % where we store most of the Data
    DisplayName     ='';
    Label           ='';
    ModificationDate=[];
    Source          = pwd;
    Tag             ='';
    Title           ='';
    User            = getenv('USER');
    UserData        =[];
 
  end % properties
  
  properties (Access=private, Hidden=true)     % internal use
    
    Private
    % Data handling: Signal, Axes, ...
    Labels= struct(); % struct: Labels.Signal, ... Labels.Axes{1:ndims}
    Axes  = {};       % cell{1:ndims} e.g. aliases
  end
  
  properties (Access=protected, Constant=true)  % shared by all instances
    Protected={'Protected','Labels','Axes','Tag','Private'} % can not be changed
  end
  
% ------------------------------------------------------------------------------
  
  methods
    function new = estruct(varargin)
    %ESTRUCT Create or convert to extended structure array.
    %  S = ESTRUCT('field1',VALUES1,'field2',VALUES2,...) creates a
    %  structure array with the specified fields and values.  The value
    %  arrays VALUES1, VALUES2, etc. must be cell arrays of the same
    %  size, scalar cells or single values.  Corresponding elements of the
    %  value arrays are placed into corresponding structure array elements.
    %  The size of the resulting structure is the same size as the value
    %  cell arrays or 1-by-1 if none of the values is a cell.
    %
    %  ESTRUCT is similar to STRUCT, but brings extended functionalities.
    %
    %  ESTRUCT(OBJ) converts the object OBJ into its equivalent
    %  structure.  The class information is lost.
    %
    %  ESTRUCT([]) creates an empty structure.
    %
    %  To create fields that contain cell arrays, place the cell arrays
    %  within a VALUE cell array.  For instance,
    %    s = estruct('strings',{{'hello','yes'}},'lengths',[5 3])
    %  creates the 1-by-1 structure
    %     s = 
    %        strings: {'hello'  'yes'}
    %        lengths: [5 3]
    %
    %  Example: s = estruct('type',{'big','little'},'color','red','x',{3 4})
    %
    %  See also isstruct, setfield, getfield, fieldnames, orderfields, 
    %  isfield, rmfield, deal, substruct, struct2cell, cell2struct.
    
      persistent id
      
      warning('off','MATLAB:structOnObject');
      new.Private.cache = []; % init cache to empty
      % handle Tag number
      if isempty(id) id=0; end
      if id > 1e6,   id=0; end % use clock
      if id <=0, 
        id = new.Date;
        id = fix(id(6)*1e4); 
      else 
        id=id+1;
      end
      new.Tag = [ 'iD' sprintf('%0.f', id) ]; % unique ID
      new.ModificationDate = new.Date;
      % add our 'static' properties so that they are equially handled by
      % subsref/subsasgn
      new.addprop('Error');           % e.g. an alias or empty
      new.addprop('Monitor');         % e.g. an alias or empty
      new.addprop('Signal');          % e.g. an alias
      new.Error = 'matlab: sqrt(this.Signal)';
      if ~nargin, return; end
      
      % collect items to store: as structures, as data files, and others
      structs = {}; % cell: will contain struct('name','value')
      
      % append arguments (not from files)
      index=1;
      while index<=nargin
        this = varargin{index}; s = [];
        if isobject(this), this = struct(this); end
        if ischar(this)
          if index<numel(varargin) && isvarname(this)   % input: name/value pair, e.g. 'par',value, ...
            s.name = this;
            s.value= varargin{index+1};
            structs{end+1} = s; this = [];
            varargin{index} = []; varargin{index+1} = []; % clear memory
            index=index+1;  % increment name/value pair
          elseif exist(this, 'file') || any(strcmp(strtok(this, ':'), {'http' 'https' 'ftp' 'file'}))
            % pass: we handle this below when assembling the object(s)
            this = [];
          else
            s.value = str2struct(this); % convert to struct ?
            s.name  = inputname(index);
            structs{end+1} = s; this = [];
            varargin{index} = [];
          end
        end % ischar
        
        if ~isempty(this)
          s.value = this;
          s.name  = inputname(index);
          structs{end+1} = s;
        end
        
        index=index+1;
      end % varargin
      
      % fill in direct name/value pairs
      for index=1:numel(structs)
        s = structs{index};
        if isempty(s.name), s.name=sprintf('%s_%i', class(s.value), index); end
        if ~isfield(new, s.name), new.addprop(s.name); end
        new.(s.name)=s.value; 
        structs{index} = []; % clear memory
      end
      numel_new = 1;
      new.Private.cache.check_requested = true; % request a check at first 'get'
      new0 = copyobj(new);
      
      % now build the final 'master' object
      for index_varg=1:numel(varargin) % loop on initial input arguments for iLoad
        arg = varargin{index_varg};
        if isempty(arg) || (~ischar(arg) && ~iscellstr(arg)), continue; end
        if ischar(arg), arg = cellstr(arg); end
        for index_arg=1:numel(arg)  % loop on input argument content when e.g. cellstr array
          this = arg{index_arg}; % must be a char or single cellstr
          if isempty(this), continue; end
          if ~isempty(dir(this)) || any(strcmp(strtok(this, ':'), {'http' 'https' 'ftp' 'file'}))
            try
              this = iLoad(this); % imported data from file goes in Data
            catch ME
              disp([ mfilename ': WARNING: failed importing ' char(this) ]);
            end
          end
          % this can now be initial char/cellstr, a struct from iLoad, or a cell of structs from iLoad
          if ~iscell(this), this = { this }; end
          for index_data=1:numel(this)
            if numel_new == 1, new.Data = this{index_data};
            else
              new1 = copyobj(new0); new1.Data = this{index_data}; this{index_data} = [];
              new = [ new new1 ]; % build array
            end % index_data 
            numel_new = numel_new+1;
          end
        end % index_arg
        
      end % index_varg
      
    end % estruct (instantiate)
    
% ------------------------------------------------------------------------------
    
    function tf = isstruct(self)
    %  ISSTRUCT True for structures.
    %      ISSTRUCT(S) returns logical true (1) if S is a structure
    %      and logical false (0) otherwise.
    %
    % See also estruct, isfield, iscell, isnumeric, isobject.
      tf = true; % isa(self, mfilename); always true
    end
    
    function tf = isfield(self, f)
    %   ISFIELD True if field is in structure array.
    %      ISFIELD(S,FIELD) returns true if the string FIELD is the name of a
    %      field in the structure array S.
    %   
    %      TF = ISFIELD(S,FIELDNAMES) returns a logical array, TF, the same size
    %      as the size of the cell array FIELDNAMES.  TF contains true for the
    %      elements of FIELDNAMES that are the names of fields in the structure
    %      array S and false otherwise.
    %   
    %      NOTE: TF is false when FIELD or FIELDNAMES are empty.
    %   
    %      Example:
    %         s = estruct('one',1,'two',2);
    %         fields = isfield(s,{'two','pi','One',3.14})
    %   
    %      See also getfield, setfield, fieldnames, orderfields, rmfield,
    %      isstruct, estruct. 
      if nargin ~= 2, return; end
      if numel(self) == 1
        tf = isfield(struct(self), f);
      else
        tf = zeros(size(self));
        for index=1:numel(self); tf(index) = isfield(self(index), f); end
      end
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
    %      See also setfield, getfield, isfield, fieldnames.
      if nargin ~= 2, return; end
      if numel(self) == 1
        if isfield(self, f)
          delete(findprop(self, f));
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
    %      See also setalias, getalias, isfield, fieldnames.
    
    % compatibility with original iData (2007-2019)
      self = rmfield(self, varargin{:});
    end 
    
    function c = struct2cell(self)
    %   STRUCT2CELL Convert structure array to cell array.
    %      C = STRUCT2CELL(S) converts the M-by-N structure S (with P fields)
    %      into a P-by-M-by-N cell array C.
    %   
    %      If S is N-D, C will have size [P SIZE(S)].
    %   
    %      Example:
    %        clear s, s.category = 'tree'; s.height = 37.4; s.name = 'birch';
    %        c = struct2cell(s); f = fieldnames(s);
    %   
    %      See also cell2struct, fieldnames.
    
      if numel(self) == 1
        c = struct2cell(struct(self));
      else
        c = arrayfun('struct2cell', self);
      end
      
    end
    
    function new = cell2struct(self, varargin)
    %   CELL2STRUCT Convert cell array to structure array.
    %      S = CELL2STRUCT(C,FIELDS,DIM) converts the cell array C into
    %      the structure S by folding the dimension DIM of C into fields of
    %      S.  SIZE(C,DIM) must match the number of field names in FIELDS.
    %      FIELDS can be a character array or a cell array of strings.
    %   
    %      Example:
    %        c = {'tree',37.4,'birch'};
    %        f = {'category','height','name'};
    %        s = cell2struct(estruct, c,f,2);
    %   
    %      See also struct2cell, fieldnames.
      if nargin < 3, d=2; end
      if nargin < 2, error([ mfilename ': cell2struct: requires at least 2 arguments' ]); end
      if numel(self) > 1
        error([ mfilename ': cell2struct(estruct, ...) only works with a single estruct.' ])
      end
      new = copyobj(self, cell2struct(varargin{:}));
    end
    
    function s = char(self)
      %  CHAR Create character array (string).
      %      S = CHAR(X) converts the structure X into a character representation
      s = class2str(self, 'eval');
    end
    
    function v = cast(self, typ)
      %   CAST  Cast object to a different data type or class.
      %     B = CAST(A,NEWCLASS) casts A to class NEWCLASS. A must be convertible to
      %     class NEWCLASS. NEWCLASS must be the name of one of the builtin data types.
      if nargin < 2, typ=''; end
      if isempty(typ), typ='double'; end
      if numel(self) == 1
        v = cast(subsref(self,struct('type','.','subs','Signal')), typ);
      else
        v = cell(size(self));
        for index=1:numel(self); v{index} = cast(self(index), typ); end
      end
    end
    
    function v = double(self)
      %   DOUBLE Convert to double precision.
      %      DOUBLE(X) returns the double precision value of the object.
      v = cast(self, 'double');
    end
    
    function v = single(self)
      %   SINGLE Convert to SINGLE precision.
      %      SINGLE(X) returns the SINGLE precision value of the object.
      v = cast(self, 'single');
    end
    
    function v = logical(self)
      %   LOGICAL Convert to logical (0=false, any other=true).
      %      LOGICAL(X) returns the SINGLE precision value of the object.
      v = cast(self, 'logical');
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
      if numel(self) > 1
        error([ mfilename ': repmat(estruct, M,N,...) only works with a single estruct.' ])
      end
      s = estruct(repmat(struct(self), varargin{:}));
    end
    
    function o = ones(self, varargin)
      % ONES   Ones array.
      %    ONES(S,N) is an N-by-N matrix of 'S' structure.
      % 
      %    ONES(S, M,N) or ONES(S, [M,N]) is an M-by-N matrix of 'S' structure.
      % 
      %    ONES(S,M,N,P,...) or ONES(S,[M N P ...]) is an M-by-N-by-P-by-... array of
      %    'S' structure.
      if nargin == 1, o=self(1); return; end
      o = repmat(self, varargin{:});
    end
    
    function z = zeros(self, varargin)
      % ZEROS  Zeros array.
      %    ZEROS(S,N) is an N-by-N matrix of empty structures.
      % 
      %    ZEROS(S,M,N) or ZEROS(S,[M,N]) is an M-by-N matrix of empty structures.
      % 
      %    ZEROS(S,M,N,P,...) or ZEROS(S,[M N P ...]) is an M-by-N-by-P-by-... array of
      %    empty structures.
      if nargin == 1, z=estruct; return; end
      z = ones(estruct, varargin{:});
    end

    function s = struct(self)
      % STRUCT Create or convert to structure array.
      %   STRUCT(OBJ) converts the object OBJ into its equivalent
      %   structure.  The class information is lost.
      if numel(self) == 1
        s = builtin('struct',self);
      else
        s = arrayfun('struct', self);
      end
    end
    
  end % methods
  
% ------------------------------------------------------------------------------
end % classdef
