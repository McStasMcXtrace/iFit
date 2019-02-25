function s = read_yaml(filename)
% data=read_yaml(filename)
%
% read_yaml Wrapper to directly read YAML files using the YAML class
%
% Input:  filename: yaml/json file (string)
% output: structure
% Example: y=read_yaml(fullfile(ifitpath, 'Data','iFit_iD80908.yaml')); isstruct(y)
%
% (c) E.Farhi, ILL. License: EUPL.
% See also: read_json, read_xml

s=[];

if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
    s.name           = 'YAML';
    s.patterns       ='';
    s.method         = mfilename;
    s.extension      = {'yaml','yml'};
    return
end

s       = YAML.read(filename);

end

