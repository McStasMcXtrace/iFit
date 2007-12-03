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

try % disable some warnings
  warn.seta = warning('off','iData:setaxis');
  warn.geta = warning('off','iData:getaxis');
  warn.get  = warning('off','iData:get');
catch
  warn = warning('off');
end

% get Signal, Error and Monitor for 'a' and 'b'
if isa(a, 'iData') & isa(b, 'iData') 
  if strcmp(op, 'combine')
    [a,b] = union(a,b); % perform combine on union
  else
    [a,b] = intersect(a,b); % perform operation on intersection
  end
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
  if p1 & ~all(m3 == 0), s3 = s3.*m3; e3=e3.*m3; end
case {'times','rdivide', 'ldivide'}
  s3 = genop(op, y1, y2); 
  m3 = genop(@times, m1, m2);
  if p1 & ~all(m3 == 0), s3 = s3.*m3; end
  e3 = sqrt(genop(@plus, (e1./s1).^2, (e2./s2).^2)).*s3;
case {'power'}
  m3 = genop(op, m1, m2); 
  s3 = genop(op, y1, y2);
  if p1 & ~all(m3 == 0), s3 = s3.*m3; end
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
if strcmp(op, 'combine')  % dimension of result might change from original object. 
                          % Can not store, thus redefine aliases as numerical values
  c = setalias(c, 'Signal', s3);
  c = setalias(c, 'Error', abs(e3));
  c = setalias(c, 'Monitor', m3);
else
  c = set(c, 'Signal', s3, 'Error', abs(e3), 'Monitor', m3);
end
c = iData_private_history(c, op, a,b);

% reset warnings
try
  warning(warn.seta);
  warning(warn.geta);
  warning(warn.get);
catch
  warning(warn);
end
