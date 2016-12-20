function endf = read_endf(filename)
% data = read_endf(filename)
%
%   import ENDF files, either using PyNE or slower/partial Matlab reader.
%
% Useful tokens for neutron scattering
% MF1/MT451 as 'info':
%   ZA    Standard material charge
%   AWR   Standard material mass 
%   AWI   Mass of the projectile in neutron units
%   EMAX  Upper limit of energy range for evaluation.
%   TEMP  Target temperature (Kelvin) for Doppler broadening.
%
% MF7/MT4 as 'thermal_inelastic':
%   LAT   Flag indicating which temperature has been used to compute a and b
%           LAT=0, the actual temperature has been used.
%           LAT=1, the constant T0 = 0.0253 eV has been used (293 K).
%   LASYM Flag indicating whether an asymmetric S(a,b) is given
%           LASYM=0, S is symmetric.
%           LASYM=1, S is asymmetric
%   B=[sigma_free, E_elastic, A=mass/mn, E_max, ~, atom_multiplicity]
%   sigma_bound = sigma_free*(A+1)^2/A^2
%
% MF7/MT2 as 'thermal_elastic':
%   SB    Bound cross section (barns) [incoherent elastic LTHR=2]
%
% Format is defined at https://www.oecd-nea.org/dbdata/data/manual-endf/endf102.pdf
% (c) E.Farhi, ILL. License: EUPL.

%
% for MF7 MT4 from ACE
% lib = pyne.ace.Library('/home/farhi/Programs/MCNP/CAB-Sab/hwtr-293/hwtr-293.ace' )
% lib.read()
% scipy.io.savemat('ace.mat',lib.tables)
% matlab> ace=load('ace.mat')
%

persistent status

% ============================ Use PyNE ? ======================================
if ~exist('status') || isempty(status) || ~isstruct(status)
  status = read_endf_pyne_present;  % private function inline below
end

if nargin == 0, return; end

if status
  endf = read_endf_pyne(filename);  % private function inline below
  if ~isempty(endf), return; end
end

% ============================ Use pure Matlab =================================
% open file and read it entirely
endf = [];

try
  content  = fileread(filename);
catch
  disp([ mfilename ': could not read file ' filename ]);
  return;
end

% open file and read it line by line
disp([ mfilename ': Opening ' filename ]);
content  = regexp(content, '\n+', 'split');
section  = [];  % current section. Starts unset.
% global data read from MF1/MT451
MF1      = [];

% read lines one by one ========================================================
for cline=content % each line is a cellstr
  if isempty(cline), continue; end
  tline = cline{1};
  if isempty(tline), continue; end
  % all ENDF lines should be 80 chars.
  if length(tline) < 80
    disp([ mfilename ': WARNING: skipping invalid ENDF line (less than 80 chars):' ])
    disp(tline)
    continue
  end
  if isempty(tline), continue; end
  tline(~isspace(tline) & ~isstrprop(tline, 'print')) = ' ';  % replace invalid chars
  MAT   = str2double(tline(67:70));
  MF    = str2double(tline(71:72));
  MT    = str2double(tline(73:75));
  NS    = str2double(tline(76:end));
  if any(isnan([MAT MF MT NS])), continue; end
  if MAT < 0 || MF < 0 || MT < 0 || NS < 0, continue; end
  
  % detect section end =========================================================
  % [MAT,MF,    0/ 0.0, 0.0, 0, 0, 0, 0] SEND
  if MT == 0 && NS == 99999
    % end of section found: store current section
    if ~isempty(section)
      % treat specific sections. Others are stored ignored.
      section0=section;
      if section.MF == 1
        endf.info = read_endf_mf1(section0);
        MF1 = endf.info;
      elseif section.MF == 7 && section.MT == 2
        % treat MF7 and add MF7/MT2 data
        if ~isfield(endf, 'thermal_elastic'), endf.thermal_elastic=struct(); end
        endf.thermal_elastic   = read_endf_mf7_mt2(section0, MF1, endf.thermal_elastic);
      elseif section.MF == 7 && section.MT == 4
        % treat MF7 and add MF7/MT4 data
        endf.thermal_inelastic = read_endf_mf7_mt4(section0, MF1);
      end
    end
    section = []; % clear memory
  end

  % read lines and store them ==================================================
  if any(MF == [1 7]) && any(MT == [0 451 2 4]) && MAT
    % first inserts spaces in between 11-char fields (FORTRAN 6F11)
    tline = [ tline(1:66) ' ' ];
    if MF~=1 % any numerical section except General Information
      i11=1:11;
      index=[ i11 67 (i11+11) 67 (i11+11*2) 67 (i11+11*3) 67 (i11+11*4) 67 67 (i11+11*5) ];
      tline = tline(index);
    end
    % WARNING: all numerics are given as [mant][+-][exp] instead of [mant]e[+-][exp]
    % e.g. 1.010000+2 -> 1.010000e+2
    tline = regexprep(tline,'([0-9.]+)([+-]{1})([0-9.]+)','$1e$2$3');
    % store the MAT_MF_MT section
  
    if MT
      if isempty(section) % new section initialised
        section.field = sprintf('MAT%i_MF%i_MT%i', MAT, MF, MT);
        section.description = read_endf_MF(MF);
        if ~isempty(section.description)
          disp([ mfilename ': File section "' section.description '" as MAT=' num2str(MAT) ' MF=' num2str(MF) ' MT=' num2str(MT) ]);
        end
        section.MF = MF;
        section.MT = MT;
        section.MAT= MAT;
        section.lines       = { tline };
      else
        section.lines{end+1} = tline;
      end
    end
  end
  
end % for cline

% ------------------------------------------------------------------------------
function h = read_endf_mf1(MF1)
  % read the MF1 MT451 General Information block and return its structure
  %
  % we treat the MF1 MT=451 case: general information.
  h = struct();
  
  if MF1.MF ~= 1,   return; end
  if MF1.MT ~= 451, return; end
  if numel(MF1.lines) < 4, return; end
  % get the 4 first lines into a matrix
  lines4 = num2cell(str2num([ MF1.lines{1:4} ]));  % 4 lines into a cell to use deal
  %  ZA,AWR,LRP,LFI,NLIB,NMOD,ELIS,STA,LIS,LISO,NFOR,AWI,EMAX,LREL,NSUB,
  %  NVER,TEMP,LDRV,NWD,NXC
  try
    [h.ZA,  h.AWR, h.LRP, h.LFI, h.NLIB,h.NMOD, ...
     h.ELIS,h.STA, h.LIS, h.LISO,d1,    h.NFOR, ...
     h.AWI, h.EMAX,h.LREL,d2,    h.NSUB,h.NVER, ...
     h.TEMP,d3,    h.LDRV,d4,    h.NWD, h.NXC]= deal(lines4{:});
  catch ME
    disp([ mfilename ': ERROR: ' MF1.description ' invalid HEAD length (MF1/MT451)' ]);
    disp(MF1.lines(1:4)')
    return  % invalid MF1 section (wrong nb of items in HEAD)
  end
  if any([d1 d2 d3 d4])
    disp([ mfilename ': WARNING: ' MF1.description ' wrong HEAD values (MF1/MT451)' ]);
    disp(MF1.lines{1:4})
  end
  l=MF1.lines{5}; 
  if length(l) < 66, h=[]; return; end  % ENDF line too short
  h.ZSYMAM = strtrim(l(1:11)); h.ALAB=strtrim(l(12:22));
  h.EDATE  = strtrim(l(23:32));h.AUTH=strtrim(l(34:66));
  l=MF1.lines{6};
  h.REF    = strtrim(l(1:22)); h.DDATE= strtrim(l(23:32)); 
  h.RDATE  = strtrim(l(34:43));h.ENDATE= strtrim(l(56:63));
  h.HSUB   = MF1.lines(7:9);
  h.COMMENTS= MF1.lines(10:end)';
  
  h.MF=MF1.MF; h.MT=MF1.MT;  h.MAT=MF1.MAT;
  h.description = MF1.description;
  h.field       = MF1.field;
  h.NSUB_string = read_endf_NSUB(h.NSUB);
  h.NLIB_string = read_endf_NLIB(h.NLIB);
  % display found item
  disp(sprintf('%s: MF=  %3i %s', mfilename, h.MF,   h.description));
  disp(sprintf('%s: NLIB=%3i %s', mfilename, h.NLIB, h.NLIB_string));
  disp(sprintf('%s: Date     %s %s %s',     mfilename, h.EDATE, h.DDATE, h.RDATE));
  disp(sprintf('%s: Material %s MAT=%i from %s at %s, release %s',  ...
    mfilename, h.ZSYMAM, h.MAT, h.AUTH, h.ALAB, h.ENDATE));
  disp(sprintf('%s: NSUB=%3i %s', mfilename, h.NSUB,h.NSUB_string ));

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

% ------------------------------------------------------------------------------
function t = read_endf_mf7_mt2(MF7, MF1, t)
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

  % we treat the MF7 MT=2 case: TSL
  % the other cases are ignored
  if MF7.MF ~= 7, return; end
  if MF7.MT ~= 2, return; end
  if numel(MF7.lines) < 2, return; end
  
  % read the header
  HEAD = num2cell(str2num([ MF7.lines{1} ]));  % HEAD line 1 into a cell to use deal
  if numel(HEAD) < 6
    disp([ mfilename ': ERROR: ' MF7.description ' invalid HEAD length (MF7/MT2 or MT4)' ]);
    disp(MF7.lines{1});
    return; 
  end % not a HEAD line
  MF7.lines(1) = []; % remove HEAD
  
  t.MAT   = MF7.MAT; t.MF=MF7.MF; t.MT=MF7.MT;
  t.field = MF7.field;
  
  if ~isempty(MF1), 
    t.ZSYMAM= MF1.ZSYMAM; t.EDATE = MF1.EDATE;
  end

  %  ENDF: [MAT, 7, 2/ ZA, AWR, LTHR, 0, 0, 0] HEAD
  [t.ZA,t.AWR,t.LTHR,d1,d2,d3]= deal(HEAD{:});
  if any([d1 d2 d3])
    disp([ mfilename ': WARNING: ' MF7.description ' wrong HEAD values (MF7/MT2)' ]);
    disp(HEAD{:})
  end
  
  % ENDF MF7 MT2 section (elastic)
  if     t.LTHR == 1  % Coherent Elastic =======================================
    t.description = [ MF7.description ' (Coherent Elastic)' ];
    % ENDF: [MAT, 7, 2/ T0,0,LT,0,NR,NP/E/S(E,T0) ] TAB1
    [ce, MF7.lines] = read_endf_TAB1(MF7.lines);
    if isempty(ce),  
      disp([ mfilename ': ERROR: ' t.description ' TAB1 aborted (MF7/MT2/LTHR1)' ]);         
      t=[]; return; 
    end % bad format
    if any(ce.head([ 2 4 ])), 
      disp([ mfilename ': WARNING: ' t.description ' wrong head values in TAB1 (MF7/MT2/LTHR1)' ]);
      disp(ce.head)
    end % bad values TAB1/HEAD
    t.LT = ce.head(3);  % nb of temperatures to read
    t.T  = ce.head(1);  % T0 (K)
    t.NP = ce.head(6);  % Number of Bragg edges given.
    if t.LT < 0 || t.T < 0 || t.NP <= 0, 
      t=[]; return; 
    end % bad format
    t.E  = ce.X(:)';    % energy axis (eV) as row
    t.S  = ce.Y(:)';    % S(E,T0) as a row
    t.INT= ce.INT(:);   % interpolation flag
    
    % ENDF: [MAT, 7, 2/ Tn,0,LI,0,NP,0/S(E,Tn) ] LIST 1:(LT+1)
    for index=1:(t.LT+1)
      if isempty(MF7.lines), break; end
      [ce, MF7.lines] = read_endf_LIST(MF7.lines);
      if isempty(ce),  
        disp([ mfilename ': ERROR: ' t.description ' LIST aborted (MF7/MT2/LTHR1)' ]);         
        t=[]; return; 
      end % bad format
      if any(ce.head([ 2 4 ])), 
        disp([ mfilename ': WARNING: ' t.description ' wrong head values in LIST (MF7/MT2/LTHR1)' ]);
        disp(ce.head)
      end % bad values LIST
      t.T  = [ t.T     ce.head(1) ];% Tn (K)
      t.INT= [ t.INT ; ce.head(3) ];% LI: interpolation flags
      t.NP = [ t.NP  ; ce.head(5) ];% NP
      t.S  = [ t.S   ; ce.B(:)' ];  % append S(E,Tn) as rows
    end
    if ~all(t.NP == t.NP(1))
      disp([ mfilename ': WARNING: NP is not constant in Section ' section.field ' ' section.description ' (MF7/MT2/LTHR1)'])
      disp(t.NP(:)')
    end
    
  elseif t.LTHR == 2  % Incoherent Elastic =====================================
    t.description = [ MF7.description ' (Incoherent Elastic)' ];
    % ENDF: [MAT, 7, 2/ T0,0,LT,0,NR,NP/Tint/W(T) ] TAB1
    [ie, MF7.lines] = read_endf_TAB1(MF7.lines);
    if isempty(ie), 
      disp([ mfilename ': ERROR: ' t.description ' TAB1 aborted (MF7/MT2/LTHR2)' ]);
      t=[]; return; 
    end % bad format
    if any(ie.head(2:4)), 
      disp([ mfilename ': WARNING: ' t.description ' wrong head values in TAB1 (MF7/MT2/LTHR2)' ]);
      disp(ie.head)
    end % bad values TAB1/HEAD
    t.SB = ie.head(1);  % bound cross section (barns)
    t.NR = ie.head(5);
    t.NP = ie.head(6);  % number of temperatures
    t.T  = ie.X(:)';    % temperature (K)
    t.W  = ie.Y(:)';    % Debye-Waller integral divided by the atomic mass (eV -1 ) 
                        %   as a function of temperature
    t.INT= ie.INT(:)';  % interpolation flag
  end
  
  % end read_endf_mf7_mt2

% ------------------------------------------------------------------------------
function t = read_endf_mf7_mt4(MF7, MF1)
  % ENDF MF7 MT4 section (inelastic)
  
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

  
  t = struct();
  % we treat the MF7 MT=4 case: TSL
  % the other cases are ignored
  if MF7.MF ~= 7, return; end
  if MF7.MT ~= 4, return; end
  if numel(MF7.lines) < 2, return; end
  
  % read header
  HEAD = num2cell(str2num([ MF7.lines{1} ]));  % HEAD line 1 into a cell to use deal
  if numel(HEAD) < 6
    disp([ mfilename ': ERROR: ' MF7.description ' invalid HEAD length (MF7/MT2 or MT4)' ]);
    disp(MF7.lines{1});
    return; 
  end % not a HEAD line
  MF7.lines(1) = []; % remove HEAD
  
  t.MAT   = MF7.MAT; t.MF=MF7.MF; t.MT=MF7.MT;
  t.field = MF7.field;

  if ~isempty(MF1)
    t.ZSYMAM=MF1.ZSYMAM; t.EDATE=MF1.EDATE;
  end
  
  t.description = [ MF7.description ' (Incoherent Inelastic)' ];
  [t.ZA,t.AWR,d1,t.LAT,t.LASYM,d2]= deal(HEAD{:});
  % LAT: Flag indicating which temperature has been used to compute a and b
  %  LAT=0, the actual temperature has been used.
  %  LAT=1, the constant T0 = 0.0253 eV has been used.
  % LASYM: Flag indicating whether an asymmetric S(a,b) is given
  %  LASYM=0, S is symmetric.
  %  LASYM=1, S is asymmetric
  if any([d1 d2])
    disp([ mfilename ': WARNING: ' MF7.description ' wrong HEAD values (MF7/MT4)' ]);
    disp(HEAD{:})
  end
  
  % [MAT,7,4/ 0,0,LLN,0,NI,NS/B(N) ] LIST
  [ii, MF7.lines] = read_endf_LIST(MF7.lines);
  if isempty(ii), 
    disp([ mfilename ': ERROR: ' t.description ' TAB1/B aborted (MF7/MT4)' ]);
    t=[]; return; 
  end % bad format
  if any(ii.head([ 1 2 4 ])), 
    disp([ mfilename ': WARNING: ' t.description ' wrong head values in LIST/B (MF7/MT4)' ]);
    disp(ii.head)
  end % bad values LIST
  t.LLN = ii.head(3); % Flag indicating the form of S(a,b) stored in the file
                      % LLN=0, S is stored directly. LLN=1, ln S is stored.
  t.NI  = ii.head(5); % Total number of items in the B(N) list. NI = 6(NS+1).
  t.NS  = ii.head(6); % Number of non-principal scattering atom types.
  t.B   = ii.B;       % List of phys. constants.
  
  % [MAT,7,4/ 0,0,0,0,NR,NB/BetaInt ] TAB2 (interpolation scheme)
  [beta, MF7.lines] = read_endf_TAB2(MF7.lines);
  if isempty(beta), 
    disp([ mfilename ': ERROR: ' t.description ' TAB2/beta aborted (MF7/MT4)' ]);
    t=[]; return; 
  end % bad format
  if any(beta.head(1:4)), 
    disp([ mfilename ': WARNING: ' t.description ' wrong head values in TAB2/beta (MF7/MT4)' ]);
    disp(beta.head)
  end % bad values TAB2
  t.NR       = beta.head(5);  % NR: Number of interpolation ranges for a particular parameter
  t.NB       = beta.head(6);  % NB: Total number of beta values given.
  t.beta_INT = beta.INT;      % Interpolation schemes
  t.beta     = zeros(1,t.NB);
  t.T        = [];
  t.LT       = t.beta;
  t.NP       = t.NB;
  t.alpha    = [];
  t.Sab      = [];  % S(a,BETAn,Tn)
  
  % read alpha,beta,T,Sab values
  for ibeta = 1:t.NB
    % [MAT,7,4/ T0,Beta0,LT,0,NR,NP/AlphaInt/S(a,BETAn,Tn) ] TAB1 Sab(T0)
    [alpha, MF7.lines] = read_endf_TAB1(MF7.lines);
    if isempty(alpha), 
      disp([ mfilename ': ERROR: ' t.description ' TAB1/alpha aborted (MF7/MT4)' ]);
      t=[]; return; 
    end % bad format
    if any(alpha.head(4)), 
      disp([ mfilename ': WARNING: ' t.description ' wrong head values in TAB1/alpha (MF7/MT4)' ]);
      disp(alpha.head)
    end % bad values TAB1
    t.T           = alpha.head(1);  % T0
    t.beta(ibeta) = alpha.head(2);  % beta(ibeta)
    t.LT(ibeta)   = alpha.head(3);  % LT: Temperature dependence flag
    t.NR(ibeta)   = alpha.head(5);  % NR
    t.NP(ibeta)   = alpha.head(6);  % NR
    t.alpha_INT(ibeta) = alpha.INT;
    t.alpha            = alpha.X(:)';          % should all be the same
    t.Sab(:, ibeta, 1) = alpha.Y; 
    % ENDF: [MAT, 7, 4/ Tn,BETAn,LT,0,NP,0/S(a,BETAn,Tn) ] LIST
    for index=2:(t.LT+1)
      [Sab, MF7.lines] = read_endf_LIST(MF7.lines);
      if isempty(Sab), 
        disp([ mfilename ': ERROR: ' t.description ' LIST/Sab aborted (MF7/MT4)' ]);
        t=[]; return; 
      end % bad format
      if any(Sab.head([ 4 6 ])), 
        disp([ mfilename ': WARNING: ' t.description ' wrong head values in LIST/Sab (MF7/MT4)' ]);
        disp(Sab.head)
      end % bad values LIST
      t.T(index)  = Sab.head(1);% Tn (K)
      t.Sab(:, ibeta, index) = Sab.B;  % append S(E,Tn)
    end
  end % beta loop (NB)
  
  % handle ln(S) storage and compute real Sab
  if t.LLN, t.Sab = exp(t.Sab); t.LLN=0; end
  
  % read Teff values
  % [MAT,7,4/ 0,0,0,0,NR,NT/Tint/Teff(T) ] TAB1
  Teff_index=6;
  if numel(t.B) >= 13 && t.B(13) == 0, Teff_index=[ Teff_index 13 ]; end
  if numel(t.B) >= 17 && t.B(17) == 0, Teff_index=[ Teff_index 17 ]; end
  t.NT = []; t.Teff_INT = []; t.Teff_T = [];
  for index=1:numel(Teff_index)
    if ~isempty(MF7.lines) && numel(t.B) >= index && t.B(index)
      [Teff, MF7.lines] = read_endf_TAB1(MF7.lines);
      if isempty(Teff), 
        disp([ mfilename ': ERROR: ' t.description ' TAB1/Teff aborted (MF7/MT4)' ]);
        t=[]; return; 
      end % bad format
      if any(Teff.head(1:4)), 
        disp([ mfilename ': WARNING: ' t.description ' wrong head values in TAB1/Teff (MF7/MT4)' ]);
        disp(Teff.head)
      end % bad values TAB1
      t.NT(index)       = Teff.head(6);
      t.Teff_INT(index) = Teff.INT;
      t.Teff_T(:,index) = Teff.X(:);  % should be t.T
      t.Teff(:,index)   = Teff.Y(:);
    end
  end

  % end read_endf_mf7_mt4


% ------------------------------------------------------------------------------
% ENDF Data Format routines to read Records (ENDF-102 2004/page 0.26-0.30)
function [TAB1, lines] = read_endf_TAB1(lines)
  % read a ENDF TAB1  record. treated lines are removed
  % includes a TAB2
  % [MAT,MF,MT/ C1, C2, L1, L2, NR, NP/xint/y(x)] TAB1
  % returns: TAB1 structure with members: head,NBT,INT,X,Y
  l0=lines;
  if isempty(lines), TAB1=[]; return; end
  [TAB1, lines] = read_endf_TAB2(lines);  % read NBT/INT with NR elements
  if isempty(TAB1), 
    disp([ mfilename ': ERROR: TAB1: wrong HEAD length' ])
    disp(lines{1})
    return; 
  end % wrong TAB1 HEAD length
  NP   = TAB1.head(6);
  % read x(N) y(N) (NP*2 items)
  n          = min(ceil(2*NP/6), numel(lines));  % nb of lines (NP/6)
  if n<1, 
    disp([ mfilename ': ERROR: TAB1: wrong content length [n NP]' ])
    disp([n NP]);
    disp(l0{1});
    disp(l0{2});
    TAB1=[]; return; 
  end
  item_lines = str2num([ lines{1:n} ]);
  NP         = min(NP*2, numel(item_lines));
  X_Y        = item_lines(1:NP);
  TAB1.X     = X_Y(1:2:end);
  TAB1.Y     = X_Y(2:2:end);
  lines(1:n) = [];

function [TAB2, lines] = read_endf_TAB2(lines)
  % read a ENDF TAB2  record. treated lines are removed
  % [MAT,MF,MT/ C1, C2, L1, L2, NR, NZ/ Z int ] TAB2
  % returns: TAB2 structure with members: head,NBT,INT
  l0=lines;
  if isempty(lines), TAB2=[]; return; end
  TAB2.head  = str2num(lines{1});  % TAB2 HEAD line into a vector
  if length(TAB2.head) < 6, 
    disp([ mfilename ': ERROR: TAB2: wrong HEAD length' ])
    disp(lines(:))
    TAB2 = [];
    return; 
  end % wrong HEAD length
  NR         = TAB2.head(5);
  lines(1)   = [];             % remove HEAD line
  % read NBT(N) INT(N) (NR*2 items)
  n          = min(ceil(NR/6), numel(lines));  % nb of lines (NR/6) remaining
  if n<1, 
    disp([ mfilename ': ERROR: TAB2: wrong content length [n NR]' ])
    disp([ n NR ]);
    disp(l0{1});
    disp(l0{2});
    TAB2=[]; return; 
  end
  item_lines = str2num([ lines{1:n} ]);
  NR         = min(NR*2, numel(item_lines));
  NBT_INT    = item_lines(1:NR);
  TAB2.NBT   = NBT_INT(1:2:end);
  TAB2.INT   = NBT_INT(2:2:end);
  lines(1:n) = [];

function [LIST, lines] = read_endf_LIST(lines)
  % read a ENDF LIST  record. treated lines are removed
  % [MAT,MF,MT/ C1, C2, L1, L2, NPL, N2/ Bn ] LIST
  % returns: LIST structure with members: head,B
  l0=lines;
  if isempty(lines), LIST=[]; return; end
  LIST.head  = str2num(lines{1});  % LIST HEAD line into a vector
  if length(LIST.head) < 6, 
    disp([ mfilename ': ERROR: LIST: wrong HEAD length' ])
    LIST=[]; return; 
  end  % wrong HEAD length
  NPL        = LIST.head(5);
  lines(1)   = [];             % remove HEAD line
  % read B(N) INT(N) (NPL items)
  n          = min(ceil(NPL/6), numel(lines));  % nb of lines (NR/6)
  if n<1, 
    disp([ mfilename ': ERROR: LIST: wrong content length [n NPL]' ])
    disp([ n NPL ]);
    disp(l0{1});
    disp(l0{2});
    LIST=[]; return; 
  end
  item_lines = str2num([ lines{1:n} ]);
  NPL        = min(NPL, numel(item_lines));
  LIST.B     = item_lines(1:NPL);
  lines(1:n) = [];

% ------------------------------------------------------------------------------
% ENDF Data Format labels
function MF_description = read_endf_MF(index)
  % get File dscription (MF label)
  persistent MF
  
  if isempty(MF)
  
    MF = cell(1,99); % allow 99 MF sections;
    MF{1}='General information';
    MF{2}='Resonance parameter data';
    MF{3}='Reaction cross sections';
    MF{4}='Angular distributions for emitted particles';
    MF{5}='Energy distributions for emitted particles';
    MF{6}='Energy-angle distributions for emitted particles';
    MF{7}='Thermal neutron scattering law data';
    MF{8}='Radioactivity and fission-product yield data';
    MF{9}='Multiplicities for radioactive nuclide production';
    MF{10}='Cross sections for radioactive nuclide production';
    MF{12}='Multiplicities for photon production';
    MF{13}='Cross sections for photon production';
    MF{14}='Angular distributions for photon production';
    MF{15}='Energy distributions for photon production';
    MF{23}='Photo- or electro-atomic interaction cross sections';
    MF{26}='Electro-atomic angle and energy distribution';
    MF{27}='Atomic form factors or scattering functions for photo-atomic interactions';
    MF{28}='Atomic relaxation data';
    MF{30}='Data covariances obtained from parameter covariances and sensitivities';
    MF{31}='Data covariances for nu(bar)';
    MF{32}='Data covariances for resonance parameters';
    MF{33}='Data covariances for reaction cross sections';
    MF{34}='Data covariances for angular distributions';
    MF{35}='Data covariances for energy distributions';
    MF{39}='Data covariances for radionuclide production yields';
    MF{40}='Data covariances for radionuclide production cross sections';
  end
  
  if index > 0 && index < numel(MF)
    MF_description = MF{index};
  else
    MF_description = '';
  end
  if isempty(MF_description), MF_description = ''; end

function NLIB_description = read_endf_NLIB(index)
  % get NLIB Library Definition (NLIB label)
  persistent NLIB
  
  if isempty(NLIB)
  
    NLIB = cell(1,50); % allow 50 MF sections;
    NLIB{1+0  } = 'ENDF/B - United States Evaluated Nuclear Data File';
    NLIB{1+1	} = 'ENDF/A - United States Evaluated Nuclear Data File';
    NLIB{1+2	} = 'JEFF - NEA Joint Evaluated File (formerly JEF)';
    NLIB{1+3	} = 'EFF - European Fusion File (now part of JEFF)';
    NLIB{1+4	} = 'ENDF/B High Energy File';
    NLIB{1+5	} = 'CENDL - China Evaluated Nuclear Data Library';
    NLIB{1+6	} = 'JENDL - Japan Evaluated Nuclear Data Library';
    NLIB{1+31	} = 'INDL/V - IAEA Evaluated Neutron Data Library';
    NLIB{1+32	} = 'INDL/A - IAEA Nuclear Data Activation Library';
    NLIB{1+33	} = 'FENDL - IAEA Fusion Evaluated Nuclear Data Library';
    NLIB{1+34	} = 'IRDF - IAEA International Reactor Dosimetry File';
    NLIB{1+35	} = 'BROND - Russian Evaluated Nuclear Data File (IAEA version)';
    NLIB{1+36	} = 'INGDB-90 - Geophysics Data';
    NLIB{1+37	} = 'FENDL/A - FENDL activation evaluations';
    NLIB{1+41	} = 'BROND - Russian Evaluated Nuclear Data File (original version)';
  end
  
  if index >= 0 && index < numel(NLIB)
    NLIB_description = NLIB{1+index};
  else
    NLIB_description = '';
  end
  if isempty(NLIB), NLIB_description = ''; end
  
function d = read_endf_NSUB(index)
  % get NSUB sub-library Definition (NSUB label)
  
  switch(index)
  case 0     ; d='Photo-Nuclear Data';
  case 1     ; d='Photo-Induced Fission Product Yields';
  case 3     ; d='Photo-Atomic Interaction Data';
  case 4     ; d='Radioactive Decay Data';
  case 5     ; d='Spontaneous Fission Product Yields';
  case 6     ; d='Atomic Relaxation Data';
  case 10    ; d='Incident-Neutron Data';
  case 11    ; d='Neutron-Induced Fission Product Yields';
  case 12    ; d='Thermal Neutron Scattering Data';
  case 113   ; d='Electro-Atomic Interaction Data';
  case 10010 ; d='Incident-Proton Data';
  case 10011 ; d='Proton-Induced Fission Product Yields';
  case 10020 ; d='Incident-Deuteron Data';
  case 20040 ; d='Incident-Alpha data';
  otherwise  ; d='';
  end

% ------------------------------------------------------------------------------
% PyNE reader routines
function status = read_endf_pyne_present
  % read_endf_pyne_present: check for availability of PyNE
  %
  % returns a flag being 1 when available.
  status = 0;

  % test for PyNE in Python
  if ismac,      precmd = 'DYLD_LIBRARY_PATH= ;';
  elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
  else precmd=''; end
  [status, result] = system([ precmd 'python -c "import pyne"' ]);
  if status ~= 0  % not present
    status = 0;
    disp([ mfilename ': warning: would make good use of PyNE.' ])
  disp('  Get it at <http://pyne.io/>.');
  disp('  Package exists for Debian/Mint/Ubuntu at <http://packages.mccode.org>.');
  else
    status = 1;
  end
  
function endf = read_endf_pyne(filename)
  % read ENDF file using PyNE
  
  if ismac,      precmd = 'DYLD_LIBRARY_PATH= ;';
  elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
  else precmd=''; end
  
  % list of python 'dict' to store in the MAT file
  dict = 'atomic_relaxation, decay, fission, info, target, projectile, resonances, thermal_elastic, thermal_inelastic, reactions';
  dict = strtrim(regexp(dict, ',', 'split')); % split as words
  % create the python script to evaluate, and a lambda function to clean keys
  tmp = tempname;
  s = { 'from pyne.endf import Evaluation ', ...
    'import scipy.io as sio ', ...
    'import re', ...
   [ 'endf=Evaluation("' filename '") ' ], ...
    'endf.read() ', ...
    'clean = lambda varStr: re.sub("^(_)","T",  re.sub("\W|^(?=\d)","_", str(varStr))  )' };
  % dict are saved to separate temporary MAT files
  % each dict key is cleaned into variable name for matlab
  for index=1:numel(dict)
    s{end+1} = sprintf([ ...
    'for key in endf.%s.keys():\n' ...
    '    c=clean(key); endf.%s[c] = endf.%s.pop(key);\n' ...
    '    if type(endf.%s[c]) == dict:\n' ...
    '        for k in endf.%s[c]: endf.%s[c][clean(k)] = endf.%s[c].pop(k)' ], ...
      dict{index}, dict{index}, dict{index}, dict{index}, dict{index}, dict{index}, dict{index});
   
    s{end+1} = sprintf('try: sio.savemat("%s_%s", endf.%s)\nexcept: None', ...
      tmp, dict{index},dict{index});
  end
  % save material and reaction list
  s{end+1} = sprintf('sio.savemat("%s_others", {"reaction_list":endf.reaction_list, "material":endf.material})', tmp);
  % write script
  fid = fopen([ tmp '.py' ],'w');
  fprintf(fid, '%s\n', s{:});
  fclose(fid);
  % now evaluate in python
  try
    [status, result] = system([ precmd 'python ' tmp '.py' ]);
  catch
    status = 127; result = 'ERROR calling Python';
  end
  disp(result);
  if status ~= 0  % error occured
    disp([ mfilename ': ERROR: could not read file ' filename ' using PyNE. Trying pure Matlab method.' ]);
    endf = [];
    return
  end
  
  % import back data from MAT file
  endf = load([ tmp '_others' ]);
  for index=1:numel(dict)
    try
      this = load([ tmp '_' dict{index} ]);
    catch ME
      disp(getReport(ME))
      this = [];
    end
    if isstruct(this) && isempty(fieldnames(this)), this = []; end
    if ~isempty(this), endf.(dict{index}) = this; end
  end
  endf.filename = filename;
  endf.pyne_script   = char(s);
  
  % delete temporary files created
  delete([ tmp '*' ]);
  
