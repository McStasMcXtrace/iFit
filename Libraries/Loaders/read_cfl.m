function s = read_cfl(filename, varargin)
% data=read_cfl(filename, options, ...) Read CFL CrysFML/FullProf crystallography file
%
% read_cfl Read CFL CrysFML/FullProf crystallography file
%
% Input:  filename: CFL CrysFML/FullProf text file (string)
% output: structure
%
% (c) E.Farhi, ILL. License: EUPL.
% See also: read_cif, read_pdb

s=[];

if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
    CFL.name            ='CFL FullProf crystallography file';
    CFL.patterns        ={'Spgr','Atom'};
    CFL.method          =mfilename;
    CFL.options         ='--headers --fortran --catenate --fast --binary --section=Atom --silent --metadata=Spgr --metadata=Cell --metadata=Spgr ';
    CFL.extension       ={'cfl'};
    CFL.postprocess     ='opencfl';
    
    s = CFL;
    return
end

% now call read_anytext with given options
if isempty(varargin)
  varargin = { '--headers --fortran --catenate --fast --binary --section=Atom --silent --metadata=Spgr --metadata=Cell --metadata=Spgr ' };
end
s       = read_anytext(filename, varargin{:});

end

