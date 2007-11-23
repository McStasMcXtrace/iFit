function c = iData_private_binary(a, b, op)
% iData_private_binary: handles binary operations
%
% Operator may apply on an iData array and:
%   a scalar
%   a vector/matrix: if its dimensionality is lower than the object, 
%     it is replicated on the missing dimensions
%   a single iData object, which is then used for each iData array element
%     operator(a(index), b)
%   an iData array, which should then have the same dimension as the other 
%     iData argument, in which case operator applies on pairs of both arguments.
%     operator(a(index), b(index))
%
% Contributed code (Matlab Central): 
%   genop: Douglas M. Schwarz, 13 March 2006

% private functions: 
%   genop: Douglas M. Schwarz, 13 March 2006
%
% for the estimate of errors, we use the Gaussian error propagation (quadrature rule), 
% or the simpler average error estimate (derivative).

% handle input iData arrays
if isa(a, 'iData') & length(a(:)) > 1
  c = [];
  if isa(b, 'iData') & length(b(:)) == length(a(:))
    for index=1:length(a(:))
      c = [ c iData_private_binary(a(index), b(index), op) ];
    end
  
  elseif isa(b, 'iData') & length(b(:)) ~= length(a(:)) & length(b(:)) ~= 1
    iData_private_error('binary',...
  [ 'Dimension of objects do not match: 1st is ' num2str(length(a(:))) ' and 2nd is ' num2str(length(b(:))) ]);
  else
    for index=1:length(a(:))
      c = [ c iData_private_binary(a(index), b, op) ];
    end
  end
  c = reshape(c, size(a));
  return
elseif isa(b, 'iData') & length(b(:)) > 1
  c = [];
  for index=1:length(b(:))
    c = [ c iData_private_binary(a, b(index), op) ];
  end
  return
end

% get Signal, Error and Monitor for 'a' and 'b'
if isa(a, 'iData') & isa(b, 'iData') 
  [a,b] = intersect(a,b); % perform operation on intersection
end
if ~isa(a, 'iData') 
  s1= a; e1=0; m1=0; p1=0;
  c = copyobj(b);
else
  s1 = get(a, 'Signal');
  e1 = get(a, 'Error');
  m1 = get(a, 'Monitor');
  c  = copyobj(a);
end
if ~isa(b, 'iData') 
  s2= b; e2=0; m2=0; p2=0;
else
  s2 = get(b, 'Signal');
  e2 = get(b, 'Error');
  m2 = get(b, 'Monitor');
end
if (all(m1==0) | all(m1==1)) & (all(m2==0) | all(m2==1)) m1=0; m2=0; end
if strcmp(op, 'combine'), p1 = 0; else p1 = 1; end

% do test on dimensionality for a vector/matrix input
% use vector duplication to fill iData dimensionality (repmat/kron)

% the 'real' signal is 'Signal'/'Monitor', so we must perform operation on this.
% then we compute the associated error, and the final monitor
% finally we multiply the result by the monitor.
if not(all(m1 == 0)) & p1, y1 = s1./m1; d1 = e1./m1; else y1=s1; d1=e1; end
if not(all(m2 == 0)) & p1, y2 = s2./m2; d2 = e2./m2; else y2=s2; d2=e2; end

% operation
switch op
case {'plus','minus','combine'}
  if strcmp(op, 'combine'), 
       s3 = genop( @plus,  s1, s2); % @plus without Monitor nomalization
  else s3 = genop( op,  y1, y2); end
  m3 = genop(@plus, m1, m2);
  e3 = sqrt(genop(@plus, d1.^2,d2.^2));
  if p1, s3 = s3.*m3; e3=e3.*m3; end
case {'times','rdivide', 'ldivide'}
  s3 = genop(op, y1, y2); 
  m3 = genop(@times, m1, m2);
  if p1, s3 = s3.*m3; end
  e3 = sqrt(genop(@plus, (e1./s1).^2, (e2./s2).^2)).*s3;
case {'power'}
  m3 = genop(op, m1, m2); 
  s3 = genop(op, y1, y2);
  if p1, s3 = s3.*m3; end
  e2logs1 = genop(@times, e2, log(s1));
  s2e1_s1 = genop(@times, s2, e1./s1);
  e3 = s3.*genop(@plus, s2e1_s1, e2logs1);
case {'lt', 'gt', 'le', 'ge', 'ne', 'eq', 'and', 'or', 'xor', 'isequal'}
  s3 = genop(op, y1, y2);
  e3 = sqrt(genop( op, d1.^2, d2.^2));
  e3 = 2*e3./genop(@plus, y1, y2); % normalize error to mean signal
  m3 = 1;
otherwise
  if isa(a,'iData'), a=a.Tag; else a=num2str(a); end
  if isa(b,'iData'), b=b.Tag; else b=num2str(b); end
  iData_private_error('binary',['Can not apply operation ' op ' on objects ' a ' and ' b '.' ]);
end

% update object
c = set(c, 'Signal', s3, 'Error', abs(e3), 'Monitor', m3);
c = iData_private_history(c, op, a,b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRIVATE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = genop(op,x,y)
%GENOP Generalized array operations.
%   GENOP(OP, X, Y) applies the function OP to the arguments X and Y where
%   singleton dimensions of X and Y have been expanded so that X and Y are
%   the same size, but this is done without actually copying any data.
%
%   OP must be a function handle to a function that computes an
%       element-by-element function of its two arguments.
%
%   X and Y can be any numeric arrays where non-singleton dimensions in one
%       must correspond to the same or unity size in the other.  In other
%       words, singleton dimensions in one can be expanded to the size of
%       the other, otherwise the size of the dimensions must match.
%
%   For example, to subtract the mean from each column, you could use
%
%       X2 = X - repmat(mean(X),size(X,1),1);
%
%   or, using GENOP,
%
%       X2 = genop(@minus,X,mean(X));
%
%   where the single row of mean(x) has been logically expanded to match
%   the number of rows in X, but without actually copying any data.
%
%   GENOP(OP) returns a function handle that can be used like above:
%
%       f = genop(@minus);
%       X2 = f(X,mean(X));

% written by Douglas M. Schwarz
% email: dmschwarz (at) urgrad (dot) rochester (dot) edu
% 13 March 2006

% This function was inspired by an idea by Urs Schwarz (no relation) and
% the idea for returning a function handle was shamelessly stolen from
% Duane Hanselman.

% Check inputs.
if ~(nargin == 1 || nargin == 3)
	error('genop:zeroInputs','1 or 3 arguments required.')
end
if ~isa(op,'function_handle')
	error('genop:incorrectOperator','Operator must be a function handle.')
end
if nargin == 1
	z = @(x,y) genop(op,x,y);
	return
end

% Compute sizes of x and y, possibly extended with ones so they match
% in length.
nd = max(ndims(x),ndims(y));
sx = size(x);
sx(end+1:nd) = 1;
sy = size(y);
sy(end+1:nd) = 1;
dz = sx ~= sy;
dims = find(dz);
num_dims = length(dims);

% Eliminate some simple cases.
if num_dims == 0 || numel(x) == 1 || numel(y) == 1
	z = op(x,y);
	return
end

% Check for dimensional compatibility of inputs, compute size and class of
% output array and allocate it.
if ~(all(sx(dz) == 1 | sy(dz) == 1))
	error('genop:argSizeError','Argument dimensions are not compatible.')
end
sz = max([sx;sy]);
z1 = op(x(1),y(1));
if islogical(z1)
	z = repmat(logical(0),sz);
else
	z = zeros(sz,class(z1));
end

% The most efficient way to compute the result seems to require that we
% loop through the unmatching dimensions (those where dz = 1), performing
% the operation and assigning to the appropriately indexed output.  Since
% we don't know in advance which or how many dimensions don't match we have
% to create the code as a string and then eval it.  To see how this works,
% uncomment the disp statement below to display the code before it is
% evaluated.  This could all be done with fixed code using subsref and
% subsasgn, but that way seems to be much slower.

% Compute code strings representing the subscripts of x, y and z.
xsub = subgen(sy ~= sz);
ysub = subgen(sx ~= sz);
zsub = subgen(dz);

% Generate the code.
indent = 2; % spaces per indent level
code_cells = cell(1,2*num_dims + 1);
for i = 1:num_dims
	code_cells{i} = sprintf('%*sfor i%d = 1:sz(%d)\n',indent*(i-1),'',...
		dims([i i]));
	code_cells{end-i+1} = sprintf('%*send\n',indent*(i-1),'');
end
code_cells{num_dims+1} = sprintf('%*sz(%s) = op(x(%s),y(%s));\n',...
	indent*num_dims,'',zsub,xsub,ysub);
code = [code_cells{:}];

% Evaluate the code.
% disp(code)
eval(code)


function sub = subgen(select_flag)
elements = {':,','i%d,'};
selected_elements = elements(select_flag + 1);
format_str = [selected_elements{:}];
sub = sprintf(format_str(1:end-1),find(select_flag));
