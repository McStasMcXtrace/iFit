classdef estruct < dynamicprops & hgsetget
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

  properties (Access=protected)
    Private
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
      warning('off','MATLAB:structOnObject');
      new.Private.cache = []; % init cache
      
      % append arguments
      index=1;
      while index<=numel(varargin)
        this = varargin{index};
        if isobject(this), this = struct(this); end
        if ischar(this) && index<numel(varargin) && isvarname(this)  % name/value pair
          if ~isfield(new, this), new.addprop(this); end
          setfield(new, this, varargin{index+1});
          index=index+1;
        elseif ischar(this) && isvarname(this) && ~isfield(new, this)
          new.addprop(this);
        elseif isstruct(this)
          % pass and treat below
        elseif ischar(this)
          s = str2struct(this); % convert to struct ?
          if isstruct(s)
            this = s; % will be used below
          else
            name=sprintf('prop_%i', index);
            if ~isfield(new, name), new.addprop(name); end
            setfield(new, name, this);
          end
        else
          name = inputname(index); 
          if isempty(name), name=sprintf('prop_%i', index); end
          if ~isfield(new, name), new.addprop(name); end
          setfield(new, name, this);
        end
        
        if isstruct(this)
          if index == 1
            new = copyobj(new, this);
          else
            new = cat(new, this);
          end
        end
        
        index=index+1;
      end
    end
    
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
    
  end
end
