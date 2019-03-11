classdef estruct < dynamicprops
%ESTRUCT Create or convert to extended structure array.
%  S = ESTRUCT('field1',VALUES1,'field2',VALUES2,...) creates a
%  structure array with the specified fields and values.  The value
%  arrays VALUES1, VALUES2, etc. must be cell arrays of the same
%  size, scalar cells or single values.  Corresponding elements of the
%  value arrays are placed into corresponding structure array elements.
%  The size of the resulting structure is the same size as the value
%  cell arrays or 1-by-1 if none of the values is a cell.
%
%  ESTRUCT is similar to STRUCT, but brings extended functionalities and can be
%  used to build sub-classes (inherit).
%
%  ESTRUCT(OBJ) converts the object OBJ into its equivalent
%  structure.  The class information is lost.
%
%  ESTRUCT([]) creates an empty structure.
%
%  ESTRUCT('filename') imports the file name and store its content into a Data property.
%
%  To create fields that contain cell arrays, place the cell arrays
%  within a VALUE cell array.  For instance,
%    s = estruct('strings',{{'hello','yes'}},'lengths',[5 3])
%  creates the 1-by-1 structure
%     s = 
%        strings: {'hello'  'yes'}
%        lengths: [5 3]
%
%  Example
%     s = estruct('type',{'big','little'},'color','red','x',{3 4})
%
%  See also isstruct, setfield, getfield, fieldnames, orderfields, 
%  isfield, rmfield, deal, substruct, struct2cell, cell2struct.

properties

    % MetaData properties
    Title
    Source = pwd;
    Creator= mfilename;
    User   = getenv('USER')
    Date   = clock;
    ModificationDate
    Command
    UserData
    Label
    DisplayName
    
    % Data properties
    Data   =[];
    Signal          % e.g. an alias
    Error           % e.g. an alias or empty
    Monitor         % e.g. an alias or empty
  end % properties
  
  properties (Access=private, Hidden=true)     % internal use
    
    Private
    % Data handling: Signal, Axes, ...
    Labels  % Labels.Signal, ... Labels.Axes{1:ndims}
    Axes    % {1:ndims} e.g. aliases
  end
  
  properties (SetAccess=private)  % can be shown, but not changed
    Tag
  end
  
  properties (Access=protected, Constant=true)  % shared by all instances
    Protected={'Labels','Axes','Tag','Private'} % can not be changed
  end
  
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
    %  Example
    %     s = estruct('type',{'big','little'},'color','red','x',{3 4})
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
      
      % append arguments
      index=1;
      while index<=numel(varargin)
        this = varargin{index};
        if isobject(this), this = struct(this); end
        if ischar(this)
          if index<numel(varargin) && isvarname(this)   % input: name/value pair, e.g. 'par',value, ...
            if ~isfield(new, this), new.addprop(this); end
            set(new, this, varargin{index+1});
            index=index+1;  % increment name/value pair
            this = [];
          elseif exist(this, 'file') || any(strcmp(strtok(this, ':'), {'http' 'https' 'ftp' 'file'}))
            try
              this1 = iLoad(this); this=[]; this.Data=this1; clear this1; % imported data from file goes in Data
            end
          else
            s = str2struct(this); % convert to struct ?
            if isstruct(s)
              this = []; this.Data = s; % will be used below, stored into Data
            end
          end
        end
        
        if isstruct(this)
          if index == 1
            new = copyobj(new, this);
          else
            new = cat(new, this);
          end
        elseif ~isempty(this)
          name = inputname(index);
          if isempty(name), name=sprintf('prop_%i', index); end
          if ~isfield(new, name), new.addprop(name); end
          set(new, name, this);
        end
        
        index=index+1;
      end
      
      new.ModificationDate = new.Date;
      
    end % estruct instantiate
    
    function tf = isstruct(self)
    %  ISSTRUCT True for structures.
    %      ISSTRUCT(S) returns logical true (1) if S is a structure
    %      and logical false (0) otherwise.
    %
    % See also estruct, isfield, iscell, isnumeric, isobject.
      tf = isa(self, mfilename);
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
      tf = isfield(struct(self), f);
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
      if nargin == 2 && isfield(self, f)
        delete(findprop(self, f));
      end
    end
    
    function [c,f] = struct2cell(self)
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
      c = struct(self);
      c = struct2cell(c);
      if nargout > 1, f = fieldnames(self); end
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
      new = copyobj(self, cell2struct(varargin{:}));
    end
    
    function [f,v] = max(self, option)
      % MAX    Largest component.
      %  MAX(X) is the biggest/largest element in X, searched recursively.
      %  [F,V] = MAX(X) returns the field name and value of the largest field.
      %
      %  [..] = MAX(X, 'numeric') returns the largest numeric field.
      %  [..] = MAX(X, 'char')    returns the largest character field.
      if nargin < 2, option = ''; end
      
      f = findfield(self, '', [ option ' biggest' ]);
      if nargout > 1
        v = get(self, f, 'link');
      end
    end
    
    function s = char(self)
      %  CHAR Create character array (string).
      %      S = CHAR(X) converts the structure X into a character representation
      s = class2str(self, 'eval');
    end
    
    function v = double(self)
      %   DOUBLE Convert to double precision.
      %      DOUBLE(X) returns the double precision value of the biggest numeric value in X.
      [~,v] = double(max(self, 'numeric'));
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
      if nargin == 1, o=self; return; end
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
    
    % pack, full, event (add listener ?), load, save, ones, zeros
    % setalias = set (does not resolve links)
    % getalias = get (does not resolve links)
  end
end
