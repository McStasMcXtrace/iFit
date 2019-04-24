function s = cell2struct(self)
% CELL2STRUCT Convert model parameters to structure
%   S = CELL2STRUCT(M) converts the model M parameters to a structure.

if numel(self) > 1
  s = {};
  for index=1:numel(self)
    s{end+1} = cell2struct(self(index));
  end
  return
end

% single model object
p = self.ParameterValues(:);
n = self.Parameters(1:numel(p));
s = cell2struct(num2cell(p),strtok(n(:)),1);