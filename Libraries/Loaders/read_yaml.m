function s = read_yaml(filename)
% data=read_yaml(filename)
%
% read_yam Wrapper to directly read YAML files using the YAML class
% Input:  filename: yaml/json file (string)
% output: structure
%
% (c) E.Farhi, ILL. License: EUPL.


s       = YAML.read(filename);

end

