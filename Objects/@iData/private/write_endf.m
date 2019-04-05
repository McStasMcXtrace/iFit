function filename = write_endf(Sab, filename)
% filename = write_endf(Sab, filename)
%
%   exports thermal scattering law sections (MF7/MT2 and 4) ENDF files. 
%   When succesful, the filename is returned, else returns empty ''.
% 
%   Currently properly exports:
% == === =============================================== ========
% MF MT  Description                                     Complete
% == === =============================================== ========
% 1  451 Descriptive data and directory                  Yes
% 7  2   Thermal elastic scattering                      Yes
% 7  4   Thermal inelastic scattering                    Yes
% == === =============================================== ========
% Other sections are ignored.
%
% Useful tokens for neutron scattering
% MF1:
%   ZA    Standard material charge
%   AWR   Standard material mass 
%   AWI   Mass of the projectile in neutron units
%   EMAX  Upper limit of energy range for evaluation.
%   TEMP  Target temperature (Kelvin) for Doppler broadening.
%
% MF7/MT4:
%   LAT   Flag indicating which temperature has been used to compute a and b
%           LAT=0, the actual temperature has been used.
%           LAT=1, the constant T0 = 0.0253 eV has been used (293 K).
%   LASYM Flag indicating whether an asymmetric S(a,b) is given
%           LASYM=0, S is symmetric.
%           LASYM=1, S is asymmetric
%   B=[sigma_free, E_elastic, A=mass/mn, E_max, ~, atom_multiplicity]
%   sigma_bound = sigma_free*(A+1)^2/A^2
%
% MF7/MT2:
%   SB    Bound cross section (barns) [incoherent elastic LTHR=2]
%
% Format is defined at https://www.oecd-nea.org/dbdata/data/manual-endf/endf102.pdf
% (c) E.Farhi, ILL. License: EUPL.

% should be used with an array of [S(alpha,beta) temperatures], or a single Sab(alpha,beta,T)

error([ mfilename ': this is not yet functional. Sorry.' ]);

% checks/assert: ZA,AWR,classical=~LASYM,temperature
%if ~isfield(Sab(1),'ZA')
%if ~isfield(endf(1),'AWR') && ~isfield(endf(1),'mass') && ~isfield(endf(1),'weight')
% must be 2D iData Sab_check would be OK ?


% initiate the findfield 'cache' by searching the mass
mass = findfield(Sab, 'mass');
if isempty(mass), 
  disp([ mfilename ': WARNING: Object ' Sab.Tag ' "' Sab.Title '" does not contain any mass information.' ])
  filename=[];
  return;
end

% ENDF compatibility flags: MF7 MT4
Sab.MAT=1; % H
Sab.MF=7;
Sab.MT=4;

Sab.ZA=101; % H in H2O
Sab.AWR=M;
Sab.LASYM=~s.classical;
Sab.Material=strtok(s.Title);
Sab.charge=Sab.ZA;
Sab.EDATE=[ 'EVAL-' upper(datestr(now,12)) ];
% ENDF compatibility flags: MF1 MT451
Sab.LRP=-1;
Sab.LFI=0;
Sab.NLIB=0;
Sab.NMOD=0;
Sab.ELIS=0;
Sab.STA=0;
Sab.LIS=0;
Sab.LISO=0;
Sab.NFOR=6;
Sab.AWI=1;
Sab.EMAX=0;
Sab.LREL=1;
Sab.NSUB=12;  % thermal
Sab.NVER=7;
Sab.TEMP=0;
Sab.LDRV=0;
Sab.NWD=53;
Sab.NXC=2;
Sab.ZSYNAM='H(H2O)';
Sab.ALAB='IKE,LANL';
Sab.AUTH='MacFarlane,Keinert,Mattes';
Sab.REF='INDC-NDS-0470';
Sab.DDATE=Sab.EDATE;
Sab.RDATE=Sab.EDATE;
Sab.ENDATE=datestr(now, 'yyyyMMDD');
Sab.HSUB={'----ENDF/B-VII.1      MATERIAL                                     ', ...
    '-----THERMAL NEUTRON SCATTERING DATA                               ', ...
    '------ENDF-6 FORMAT                                                '};
Sab.COMMENTS={''};

% search for MF1/MT451 General Information

% search for MF7/MT2 temperatures

% search for MF7/MT4 temperatures

% open the file. If exists, return (no overwrite).

% write MF1/MT451

% write MF7/MT2

% write MF7/MT4

% ------------------------------------------------------------------------------

function info = search_endf_info()  
% search for MF1/MT451 items and constructs an 'info' structure

  % PyNE              ENDF
  % =====================================
  % description       COMMENT
  % reference         REF
  % format            NFOR
  % date_distribution DDATE
  % date_release      RDATE
  % author            AUTH
  % derived           LDRV
  % sublibrary        NSUB  (NSUB_string)
  % library           NLIB  (NLIB_string)
  % energy_max        EMAX
  % date              EDATE
  % modification      NMOD
  % laboratory        ALAB
  % identifier        HSUB
  % date_entry        ENDATE

function thermal_elastic = search_endf_thermal_elastic()
% search for MF7/MT2 items and constructs a 'thermal_elastic' structure

  % ENDF MF7 MT2 section (elastic)
  %   LTHR=1  coherent
  %   LTHR=2  incoherent

  % PyNE              ENDF
  % =====================================
  % S.<T>               S
  %      .x             E
  %      .y             S
  %      .n_pairs       NP
  %      .interp        INT
  % type                LTHR=1: coherent; LTHR=2: incoherent

function thermal_inelastic = search_endf_thermal_inelastic()
% search for MF7/MT4 items and constructs a 'thermal_inelastic' structure

  % PyNE              ENDF
  % =====================================
  % B                 B
  % teff              Teff
  % temperature       T
  % scattering_law    Sab
  % ln_S_             LLN
  % beta              beta
  % symmetric         ~LASYM
  % temperature_used  LAT=0: actual temperature; LAT=1: 0.0253 eV (i.e. 293K)
  % alpha             alpha
  % num_non_principal NS

% ------------------------------------------------------------------------------
% ENDF Data Format routines to write Records (ENDF-102 2004/page 0.26-0.30)
function NS = write_endf_TAB1(fid, TAB1, MAT, MF, MT, NS)
  % write a ENDF TAB1  record.
  % includes a TAB2
  % [MAT,MF,MT/ C1, C2, L1, L2, NR, NP/xint/y(x)] TAB1
  % input:
  %   TAB1: a structure with members: head,NBT,INT,X,Y
  % returns: 
  %   NS: new line index or -1 (error)
  
  if numel(TAB1.head) < 6, NS=-1; return; end
  NP = numel(TAB1.X);
  TAB1.head(6) = NP;
  NS = write_endf_TAB2(fid, TAB1, MAT, MF, MT, NS);  % write NBT/INT with NR elements
  if NS < 0, return; end
  X_Y  = zeros(1,2*NP);
  if numel(TAB1.Y) ~= NP, NS=-1; return; end
  X_Y(1:2:end) = TAB1.X;
  X_Y(2:2:end) = TAB1.Y;
  NS = write_endf_vector(fid, X_Y, MAT, MF, MT, NS);

function NS = write_endf_TAB2(fid, TAB2, MAT, MF, MT, NS)
  % write a ENDF TAB2  record. 
  % [MAT,MF,MT/ C1, C2, L1, L2, NR, NZ/ Z int ] TAB2
  % input:
  %   TAB2: a structure with members: head,NBT,INT
  % returns: 
  %   NS: new line index or -1 (error)
  
  if numel(TAB2.head) < 6, NS=-1; return; end
  NR = numel(TAB2.NBT);
  TAB2.head(5) = NR;
  if numel(TAB2.INT) ~= NR, NS=-1; return; end
  % NZ=TAB2.head(6) can be the number of further LIST calls (NZ)
  % assemble NBT_INT
  NBT_INT=zeros(1,2*NZ);
  NBT_INT(1:2:end) = TAB2.NBT;
  NBT_INT(2:2:end) = TAB2.INT;
  TAB2.B = NBT_INT;
  % write TAB2 as a LIST
  NS = write_endf_LIST(fid, TAB2, MAT, MF, MT, NS);
  % and further calls can write LIST 1:(TAB2.head(6)). 

function NS = write_endf_LIST(fid, LIST, MAT, MF, MT, NS)
  % write a ENDF LIST  record.
  % [MAT,MF,MT/ C1, C2, L1, L2, NPL, N2/ Bn ] LIST
  % input:
  %   LIST: a structure with members: head,B
  % returns: 
  %   NS: new line index or -1 (error)
  
  if numel(LIST.head) < 6, NS=-1; return; end
  % compute how many lines should be written
  NPL= numel(LIST.B);
  NL = floor(NPL/6.0);
  % LIST.head(5) = NPL; commented-out so that it can be used as a TAB2 with length 2*NPL
  % write header line [ C1, C2, L1, L2, NPL, N2 ]
  % HEAD: C1=ZA, C2=AWR, NPL=numel(B)
  fprintf(fid, '%11g%11g%11i%11i%11i%11i', LIST.head(1:6));
  fprintf(fid, '%4i%2i%3i%5i\n', MAT,MF,MT,NS);
  NS=NS+1; if NS > 99999, NS=0; end
  NS = write_endf_vector(fid, LIST.B, MAT, MF, MT, NS)
  
function NS = write_endf_vector(fid, VECT, MAT, MF, MT, NS)
  % write a ENDF single vector.
  % input:
  %   VECT: a vector
  % returns: 
  %   NS: new line index or -1 (error)
  Bindex0=1;
  % write the LIST
  for index=1:NL
    Bindex1 = Bindex0+6;
    fprintf(fid, '%11g', VECT(Bindex0:min(NPL,Bindex1)));
    if NPL<Bindex1
      % fill incomplete line with spaces
      fprintf(fid, '%c', ones(Bindex1-NPL,11)*' ');
    end
    fprintf(fid, '%4i%2i%3i%5i\n', MAT,MF,MT,NS);
    NS=NS+1; if NS > 99999, NS=0; end
    Bindex0=Bindex1;
  end
  


