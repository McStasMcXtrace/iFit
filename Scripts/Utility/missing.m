function s = missing(obj1, obj2)
% MISSING Compare class methods.
%   MISSING(obj1, obj2) compare methods from two classes and show differences.

  s = [];
  if ~nargin, return; end

  if ischar(obj1), obj1 = feval(obj1); end

  if nargin == 1
    s=methods(obj1);
  end
  
  if ischar(obj2), obj2 = feval(obj2); end
  
  e=methods(obj1);
  d=methods(obj2);

  disp([ mfilename ': Comparing classes ' class(obj1) ' and ' class(obj2) '.' ]);
  disp('--------------------------------------------------------------------------------')
  
  disp([ mfilename ': Methods from ' class(obj1) ' missing in ' class(obj2) ])
  s.([ 'missing_in_' class(obj2) ]) = missing_methods(e, d);
  
  disp(' ')
  disp([ mfilename ': Methods from ' class(obj2) ' missing in ' class(obj1) ])
  s.([ 'missing_in_' class(obj1) ]) = missing_methods(d, e);


% ----------------------
function s = missing_methods(e, d)

  s = {};
  for m=e'
    if ~any(strcmp(m{1}, d))
      s{end+1} = m{1};
      disp([ '  ' m{1} ])
    end
  end
