function s = cell2struct(self, p, varargin)
% CELL2STRUCT Convert model parameters to structure
%   S = CELL2STRUCT(M) converts the model M parameters to a structure.
%   The current parameters are used, or guessed when not set yet.
%
%   S = CELL2STRUCT(M,P) converts the model M parameters to a structure,
%   using specified parameter values P, which should be given as a vector
%   (with values matching parameter names order) or a string such as
%   'guess' or current'.
%
%   S = CELL2STRUCT(M,'guess') guess parameter values, and then converts
%   the model M parameters to a structure.

if nargin < 2, p = []; end
if numel(self) > 1
  s = {};
  for index=1:numel(self)
    s{end+1} = cell2struct(self(index), p, varargin{:});
  end
  return
end

% single model object
if isempty(p)
  p = self.ParameterValues(:);
elseif ischar(p) && any(strcmp(p, {'guess','current'}))
  p = feval(self, p, varargin{:});
end
n = self.Parameters(1:numel(p));
s = cell2struct(num2cell(p),strtok(n(:)),1);