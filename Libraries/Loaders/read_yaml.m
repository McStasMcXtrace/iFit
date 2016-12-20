function s = read_yaml(filename)
% data=read_yaml(filename)
%
% read_yaml Wrapper to directly read YAML files using the YAML class
% Input:  filename: yaml/json file (string)
% output: structure
%
% (c) E.Farhi, ILL. License: EUPL.
% See also: read_json

s=[];
if nargin == 0, return; end

s       = YAML.read(filename);

end

