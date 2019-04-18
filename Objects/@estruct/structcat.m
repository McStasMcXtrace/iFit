function Res = structcat(A,varargin)
% STRUCTCAT concatenate (merge) structures
%
% Res=structcat(A,B)
%   Recursively merges fields and subfields of structures A and B to result structure Res
%   Simple recursive algorithm merges fields and subfields of two structures
%
% Res=structcat(A,B,'or')
%   The result has fields which are either in A and B (union). This is the default.
% Res=structcat(A,B,'and')
%   The result has fields which are both in A and B (intersection)
% Res=structcat(A,B,'xor')
%   The result has fields which are either in A and B but not both (difference)
%
% Example:
%   A.field1=1;
%   A.field2.subfield1=1;
%   A.field2.subfield2=2;
% 
%   B.field1=1;
%   B.field2.subfield1=10;
%   B.field2.subfield3=30;
%   B.field3.subfield1=1;
% 
%   C=structcat(A,B);
%
% Version: $Date$ $Version$ $Author$
% Contribution from Igor Kaufman, 02 Dec 2011, BSD
% <http://www.mathworks.com/matlabcentral/fileexchange/34054-merge-structures>

op='or'; % default operator

for index=1:numel(varargin)
  if ischar(varargin{index}) || isa(varargin{index}, 'function_handle')
    op=varargin{index}; 
    varargin(index)=[]; break;
  end
end

Res=[]; 
if nargin>0
    Res=A;  % we start with initial structure
end
if nargin ==1, return; end

if numel(varargin) > 1
  for index = 1:numel(varargin)
    Res = structcat(A, varargin{index}, op);
  end
  return
else B = varargin{1};
end

if nargin==1 || isstruct(B)==0
    return;
end

fna= fieldnames(A);
fnb= fieldnames(B);
fn = unique([fna ; fnb ]);

for i=1:length(fn) % loop on B fields
  s=char(fn(i));
  isfieldA = isfield(A,s);
  isfieldB = isfield(B,s);
  isfieldR = isfield(Res,s);
  fieldA=[]; fieldB=[]; fieldR=[];
  if isfieldA
    fieldA=getfield(A,s);
  end
  if isfieldB
    fieldB=getfield(B,s);
  end
  if isfieldR
    fieldR=getfield(Res,s);
  end

  try
    addme = feval(op, isfieldA, isfieldB); % apply operator 'op'='or','and','xor'
  catch
    addme = false;
  end
  if addme
    if isstruct(fieldA) && isstruct(fieldB)
      Res=setfield(Res,s,structcat(fieldA, fieldB, op));  % recursive inside structures
    elseif ~isfieldR
      Res.addprop(s);
      if isempty(fieldA) 
        Res=setfield(Res,s,fieldB);
      elseif isempty(fieldA) 
        Res=setfield(Res,s,fieldA);
      else % in A and B: catenate
        if isnumeric(fieldA) && isnumeric(fieldB) && isvector(fieldA) && isvector(fieldB)
          Res=setfield(Res,s,[ fieldA(:); fieldB(:) ]);
        else
          Res=setfield(Res,s,{ fieldA, fieldB });
        end
      end
    elseif isempty(fieldA)
      Res=setfield(Res,s,fieldB);
    elseif isempty(fieldB)
      Res=setfield(Res,s,fieldA);
    elseif ~isempty(fieldB) % already in Res (from A)
      if isnumeric(fieldR) && isnumeric(fieldB) && isvector(fieldR) && isvector(fieldB)
        Res=setfield(Res,s,[ fieldR(:); fieldB(:) ]);
      else
        Res=setfield(Res,s,{ fieldR, fieldB });
      end
    end
  elseif isfield(Res,s)
    Res=rmfield(Res,s);
  end

end
history(A, mfilename, A, varargin{:});

% union:     OR:  isfield(A,s) || isfield(B,s)
% intersect: AND: isfield(A,s) && isfield(B,s)
% diff:      XOR: (isfield(A,s) && ~isfield(B,s)) || (~isfield(A,s) && isfield(B,s))