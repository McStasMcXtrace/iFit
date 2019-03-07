function s = read_hzb_tas(filename, varargin)
% data=read_hzb_tas(filename, options, ...) Read HZB FLEXX FlafCone data
%
% read_hzb_tas read HZB FLEXX FlafCone text file format
%
% Input:  filename: HZB FLEXX FlafCone Data text file (string)
% output: structure
%
% (c) E.Farhi, ILL. License: EUPL.
% See also: read_llb_tas, read_anytext, read_spec

s=[];

if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
    HZB_FELXX_Flat.name    ='HZB FLEXX FlafCone';
    HZB_FELXX_Flat.patterns={'RRRR','AAAA','VVVV','DATA_:','flat:','POSQE:','PARAM:'};
    HZB_FELXX_Flat.options =[ '--headers --fortran --catenate --fast --binary --silent ' ...
                        '--section=PARAM --section=VARIA --section=ZEROS --section=DATA ' ...
                        '--section=STEPS --metadata=KFIX ' ...
                        '--metadata=LOCAL --metadata=USER --metadata=EXPNO --metadata=DATE ' ...
                        '--metadata=INSTR --metadata=COMND --metadata=TITLE --metadata=flat --metadata=PNT'];
    HZB_FELXX_Flat.method  =mfilename;
    
    s = HZB_FELXX_Flat;
    return
end

% now call read_anytext with given options

if isempty(varargin)
  varargin = { '--headers --fortran --catenate --fast --binary --makerows=FFFF --makerows=JJJJ --makerows=IIII --silent ' };
end
s       = read_anytext(filename, varargin{:});

end

