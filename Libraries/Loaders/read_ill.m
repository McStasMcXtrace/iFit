function s = read_ill(filename, varargin)
% data=read_ill(filename, options, ...) Read ILL data
%
% read_ill read Institut Laue-Langevin text file format
%
% Input:  filename: ILL Data text file (string)
% output: structure
% Example: y=read_ill(fullfile(ifitpath, 'Data','ILL_IN6.dat')); isstruct(y)
%
% $
% See also: read_llb_tas, read_anytext, read_spec

s=[];

if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
    ILL_normal.name       ='ILL Data (normal integers)';
    ILL_normal.patterns   ={'RRRR','AAAA','FFFF','IIII'};
    ILL_normal.options    ='--headers --fortran --catenate --fast --binary --makerows=IIII --makerows=FFFF --silent ';
    ILL_normal.method     =mfilename;
    
    ILL_integers.name       ='ILL Data (large integers)';
    ILL_integers.patterns   ={'RRRR','AAAA','FFFF','JJJJ'};
    ILL_integers.options    ='--headers --fortran --catenate --fast --binary --makerows=JJJJ --makerows=FFFF --silent ';
    ILL_integers.method     =mfilename;
    
    ILL_float.name          ='ILL Data (floats only)';
    ILL_float.patterns      ={'RRRR','AAAA','FFFF'};
    ILL_float.options       ='--headers --fortran --catenate --fast --binary --makerows=FFFF --silent ';
    ILL_float.method        =mfilename;
    
    ILL_general.name       ='ILL Data (general)';
    ILL_general.patterns   ={'RRRR','AAAA','SSSS'};
    ILL_general.options    ='--headers --fortran --catenate --fast --binary --makerows=FFFF --makerows=JJJJ --makerows=IIII --silent ';
    ILL_general.method     =mfilename;
    
    ILL_TAS_pol.name       ='ILL TAS Data (polarized)';
    ILL_TAS_pol.patterns   ={'PAL','POSQE:','PARAM:','DATA_:','USER_:'};
    ILL_TAS_pol.options    =['--fast --binary --headers --silent --fortran=0 --catenate ' ...
                        '--section=PARAM --section=VARIA --section=ZEROS --section=DATA ' ...
                        '--section=POLAN --section=STEPS ' ...
                        '--metadata=LOCAL --metadata=USER --metadata=EXPNO --metadata=DATE ' ...
                        '--metadata=INSTR --metadata=COMND --metadata=TITLE --metadata=MULTI --metadata=PNT'];
    ILL_TAS_pol.method     =mfilename;
    ILL_TAS_pol.postprocess='load_ill_tas'; % load_ill_tas
    
    ILL_TAS.name       ='ILL TAS Data';
    ILL_TAS.patterns   ={'POSQE:','PARAM:','DATA_:','USER_:'};
    ILL_TAS.options    =['--fast --binary --headers --silent --fortran=0 --catenate ' ...
                        '--section=PARAM --section=VARIA --section=ZEROS --section=DATA ' ...
                        '--section=STEPS ' ...
                        '--metadata=LOCAL --metadata=USER --metadata=EXPNO --metadata=DATE --metadata=DATA ' ...
                        '--metadata=INSTR --metadata=COMND --metadata=TITLE --metadata=MULTI --metadata=PNT '];
    ILL_TAS.method      =mfilename;
    ILL_TAS.postprocess ='load_ill_tas'; % load_ill_tas
    
    s = { ILL_normal, ILL_integers, ILL_float, ILL_general, ILL_TAS_pol, ILL_TAS };
    return
end

% now call read_anytext with given options

if isempty(varargin)
  varargin = { '--headers --fortran --catenate --fast --binary --makerows=FFFF --makerows=JJJJ --makerows=IIII --silent ' };
end
s       = read_anytext(filename, varargin{:});

end

