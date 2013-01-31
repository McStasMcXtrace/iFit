function a = mlock(a, varargin)
% b = mlock(s, parameters, ...) : sets parameter lock (fix) for further fit using the model 's'
%
%   @iFunc/mlock lock model parameters during further fit process.
%     to lock/fix a parameter model, use mlock(model, parameter)
%
%   mlock(model, {Parameter1, Parameter2, ...})
%     lock/fix parameter for further fits
%   mlock(model)
%     display fixed parameters
%
% input:  s: object or array (iFunc)
% output: b: object or array (iFunc)
% ex:     b=mlock(a,'Intensity');
%
% Version: $Revision: 1.1 $
% See also iFunc, iFunc/fits, iFunc/munlock, iFunc/xlim

% calls subsasgn with 'fix' for each parameter given

% handle array of objects
if numel(a) > 1
  for index=1:numel(a)
    a(index) = feval(mfilename, a(index), varargin{:});
  end
  if nargout == 0 && ~isempty(inputname(1)) % update array inplace
    assignin('caller', inputname(1), a);
  end
  return
end

if nargin == 1 % display fixed parameters
  count = 0;
  if ~isempty(a.Parameters)
    for p=1:length(a.Parameters)
      [name, R] = strtok(a.Parameters{p}); % make sure we only get the first word (not following comments)
      
      if (length(a.Constraint.fixed) >=p && a.Constraint.fixed(p)) || ...
         (length(a.Constraint.set) >=p && ~isempty(a.Constraint.set{p})) 
        if (count == 0)
          fprintf(1, 'Locked/Fixed Parameters in %s %iD model "%s"\n', a.Tag, a.Dimension, a.Name)
        end
        count = count +1;
        if length(a.Constraint.set) >=p && ~isempty(a.Constraint.set{p}) && ...
          (ischar(a.Constraint.set{p}) || isa(a.Constraint.set{p}, 'function_handle'))
          name = [ name ' set from ' char(a.Constraint.set(p)) ];
        end
        fprintf(1,'%20s (fixed)\n', name);
      end
    end
  end
  if count == 0
    fprintf(1, 'No locked/fixed Parameters in %s %iD model "%s"\n', a.Tag, a.Dimension, a.Name);
  end
  return
end

% handle multiple parameter name arguments
if length(varargin) > 1
  for index=1:length(varargin)
    a = feval(mfilename, a, varargin{index});
  end
  return
else
  name = varargin{1};
end

% now with a single input argument
if ~ischar(name) && ~iscellstr(name)
  error([ mfilename ': can not lock model parameters with a Parameter name if class ' class(name) ' in iFunc model ' a.Tag '.' ]);
end

if ischar(name), name = cellstr(name); end
% now with a single cellstr
for index=1:length(name)
  s = struct('type', '.', 'subs', name{index});
  a = subsasgn(a, s, 'fix');
end

if nargout == 0 && ~isempty(inputname(1))
  assignin('caller',inputname(1),a);
end

