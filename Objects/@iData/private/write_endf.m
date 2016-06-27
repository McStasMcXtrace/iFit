function filename = write_endf(endf, filename)
% filename = write_endf(endf, filename)
%
%   exports thermal scattering law sections (MF7/MT2 and 4) ENDF files. 
%   This writer is rather slow, but makes a full check of the ENDF file.
%   When succesful, the filename is rturned, else returns empty ''.
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

% open the file. If exists, return (no overwrite).

% get data which can be defined in non EBDF frame
% mass, temperatures, 'charge', symmetric, 

% checks/assert: ZA,AWR,classical=~LASYM,temperature
if ~isfield(endf(1),'ZA')
if ~isfield(endf(1),'AWR') && ~isfield(endf(1),'mass') && ~isfield(endf(1),'weight')

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

% write MF1/MT451

% write MF7/MT2

% write MF7/MT4

% ------------------------------------------------------------------------------
















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
sections = {};  % all found sections
section  = [];  % current section. Starts unset.
% global data read from MF1/MT451
ZSYNAM   = [];
EDATE    = [];
AWR      = [];
MF1      = [];

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
  % detect section end
  % [MAT,MF,    0/ 0.0, 0.0, 0, 0, 0, 0] SEND
  if MT == 0 && NS == 99999
    % end of section found: store current section
    if ~isempty(section)
      % treat specific sections. Others are stored as raw import data (char lines)
      section0=section;
      if section.MF == 1
        section = write_endf_mf1(section0);
        if ~isempty(section) && isfield(section,'EDATE'),  EDATE  = section(1).EDATE; end
        if ~isempty(section) && isfield(section,'ZSYNAM'), ZSYNAM = section(1).ZSYNAM; end
        if ~isempty(section) && isfield(section,'AWR'),    AWR    = section(1).AWR; end
        MF1 = section(1);
      elseif section.MF == 7
        % treat MF7 and add MF1/MT451 data
        section = write_endf_mf7(section0, ZSYNAM, EDATE, MF1);
        
      end
      % treat specific sections -> FAILED. we restore raw, and display message
      if isempty(section)
        section = section0;
        disp([ mfilename ': Failed interpreting Section ' section.field ' ' section.description '. Storing as raw char.' ])
      end
      sections{end+1} = section;
    end
    section = [];
  end
  
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
  
  if MF && MT && MAT
    if isempty(section) && MT
      section.field = sprintf('MAT%i_MF%i_MT%i', MAT, MF, MT);
      section.description = write_endf_MF(MF);
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

endf = sections;

% ------------------------------------------------------------------------------
function h = write_endf_mf1(MF1)
  % read the MF1 MT451 General Information block and return its structure
  %
  % we treat the MF1 MT=451 case: general information.
  % the cases MF=1 MT=452,455,456, 458 (fissions) are ignored
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
  h.ZSYNAM = strtrim(l(1:11)); h.ALAB=strtrim(l(12:22));
  h.EDATE  = strtrim(l(23:32));h.AUTH=strtrim(l(34:66));
  l=MF1.lines{6};
  h.REF    = strtrim(l(1:22)); h.DDATE= strtrim(l(23:32)); 
  h.RDATE  = strtrim(l(34:43));h.ENDATE= strtrim(l(56:63));
  h.HSUB   = MF1.lines(7:9);
  h.COMMENTS= MF1.lines(10:end)';
  
  h.MF=MF1.MF; h.MT=MF1.MT;  h.MAT=MF1.MAT;
  h.description = MF1.description;
  h.field       = MF1.field;
  h.NSUB_string = write_endf_NSUB(h.NSUB);
  h.NLIB_string = write_endf_NLIB(h.NLIB);
  % display found item
  disp(sprintf('%s: MF=  %3i %s', mfilename, h.MF,   h.description));
  disp(sprintf('%s: NLIB=%3i %s', mfilename, h.NLIB, h.NLIB_string));
  disp(sprintf('%s: Date     %s %s %s',     mfilename, h.EDATE, h.DDATE, h.RDATE));
  disp(sprintf('%s: Material %s MAT=%i from %s at %s, release %s',  ...
    mfilename, h.ZSYNAM, h.MAT, h.AUTH, h.ALAB, h.ENDATE));
  disp(sprintf('%s: NSUB=%3i %s', mfilename, h.NSUB,h.NSUB_string ));

% ------------------------------------------------------------------------------
function t    = write_endf_mf7(MF7, ZSYNAM, EDATE, MF1)
  % read the MF7 MT2 and MT4 "TSL" sections and return its structure
  %
  % FILE 7. THERMAL NEUTRON SCATTERING LAW DATA
  %  File 7 contains neutron scattering data for the thermal neutron energy range
  %  (E<5 eV) for moderating materials. Sections are provided for elastic (MT=2) and
  %  inelastic (MT=4) scattering. Starting with ENDF/B-VI, File 7 is complete in
  %  itself, and Files 3 and 4 are no longer required to obtain the total scattering
  %  cross section in the thermal energy range.
  
  % uses:
  % TAB1 -> head,NBT,INT,X,Y (includes TAB2)
  % TAB2 -> head,NBT,INT
  % LIST -> head,B

  t = struct();
  
  % we treat the MF7 MT=4 case: TSL
  % the other cases are ignored
  if MF7.MF ~= 7, return; end
  if MF7.MT ~= 2 && MF7.MT ~= 4, return; end
  if numel(MF7.lines) < 2, return; end

  HEAD = num2cell(str2num([ MF7.lines{1} ]));  % HEAD line 1 into a cell to use deal
  if numel(HEAD) < 6
    disp([ mfilename ': ERROR: ' MF7.description ' invalid HEAD length (MF7/MT2 or MT4)' ]);
    disp(MF7.lines{1});
    return; 
  end % not a HEAD line
  MF7.lines(1) = []; % remove HEAD
  
  t.MAT   = MF7.MAT; t.MF=MF7.MF; t.MT=MF7.MT;
  t.field = MF7.field;
  t.ZSYNAM=ZSYNAM; t.EDATE=EDATE;
  if ~isempty(MF1), t.DescriptiveData = MF1; end
  if MF7.MT == 2 % Incoherent/Coherent Elastic Scattering
    %  ENDF: [MAT, 7, 2/ ZA, AWR, LTHR, 0, 0, 0] HEAD
    [t.ZA,t.AWR,t.LTHR,d1,d2,d3]= deal(HEAD{:});
    if any([d1 d2 d3])
      disp([ mfilename ': WARNING: ' MF7.description ' wrong HEAD values (MF7/MT2)' ]);
      disp(HEAD{:})
    end
    
    [MF7,t] = write_endf_mf7_mt2(MF7, t);
    
  elseif MF7.MT == 4 % Incoherent Inelastic Scattering
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

    [MF7,t] = write_endf_mf7_mt4(MF7, t);
    
  end
  t=write_endf_mf7_array(t);
  
  % end write_endf_mf7

function t0=write_endf_mf7_array(t)
  t0 = [];
  disp(sprintf('%s: MF=%3i   MT=%i TSL T=%s', mfilename, t.MF, t.MT, mat2str(t.T)));
  for index=1:numel(t.T)
    t1     = t;
    t1.T   = t.T(index);
    t1.Title = [ t1.ZSYNAM ' T=' num2str(t1.T) ' [K] ' t1.description];
    if t1.MT == 2 % Incoherent/Coherent Elastic Scattering
      t1.NP  = t.NP(index);
      t1.S   = t.S(index,:);
      t1.INT = t.INT(index);
      t1.Label = [ t1.ZSYNAM ' T=' num2str(t1.T) ' [K] TSL elastic' ];
    elseif t1.MT == 4 % Incoherent Inelastic Scattering
      t1.Sab = t.Sab(:,:,index);
      t1.Teff     = t.Teff(index,:);
      t1.Label = [ t1.ZSYNAM ' T=' num2str(t1.T) ' [K] TSL inelastic' ];
    end
    t0 = [ t0 t1 ];
  end
% ------------------------------------------------------------------------------
function [MF7,t] = write_endf_mf7_mt2(MF7, t)
  % ENDF MF7 MT2 section (elastic)
  if     t.LTHR == 1  % Coherent Elastic
    t.description = [ MF7.description ' (Coherent Elastic)' ];
    % ENDF: [MAT, 7, 2/ T0,0,LT,0,NR,NP/E/S(E,T0) ] TAB1
    [ce, MF7.lines] = write_endf_TAB1(MF7.lines);
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
      [ce, MF7.lines] = write_endf_LIST(MF7.lines);
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
    
  elseif t.LTHR == 2  % Incoherent Elastic
    t.description = [ MF7.description ' (Incoherent Elastic)' ];
    % ENDF: [MAT, 7, 2/ T0,0,LT,0,NR,NP/Tint/W(T) ] TAB1
    [ie, MF7.lines] = write_endf_TAB1(MF7.lines);
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
  
  % end write_endf_mf7_mt2

% ------------------------------------------------------------------------------
function NS = write_endf_mf7_mt4(fid, t, MAT, MF, MT, NS)
  % ENDF MF7 MT4 section (inelastic)
  
  % check incoming data and set the unset fields
  if ~isfield(t,'B')   return; end
  if ~isfield(t,'Sab') return; end
  
  % assemble the MF7/MT4 header and write it
  % LAT: Flag indicating which temperature has been used to compute a and b
  % LASYM: Flag indicating whether an asymmetric S(a,b) is given
  % ZA,AWR: Standard material charge and mass parameters
  if ~isfield(t,'ZA')    t.ZA=1; end
  if ~isfield(t,'AWR')   t.AWR=1; end
  if ~isfield(t,'LAT')   t.LAT=0; end     % temperature must have been defined
  if ~isfield(t,'LASYM') t.LASYM=0; end   % symmetric/classical
  VECT = [ t.ZA t.AWR 0 t.LAT t.LASYM 0];
  NS = write_endf_vector(fid, VECT, MAT, MF, MT, NS);
  % assemble the LIST header and write it
  % LLN: Flag indicating the form of S(a,b) stored in the file
  % NI:  Total number of items in the B(N) list. NI = 6(NS+1).
  t.NI = numel(t.B);
  t.NS = t.NI/6 - 1;  % Number of non-principal scattering atom types.
  if ~isfield(t,'LLN') t.LLN=0; end 
  t.head = [ 0 0 t.LLN 0 t.NI t.NS ];
  % [MAT,7,4/ 0,0,LLN,0,NI,NS/B(N) ] LIST
  NS = write_endf_LIST(fid, t, MAT, MF, MT, NS);
  
  % now write the alpha, beta and Sab(beta,alpha,T)
  
  if ~isfield(t,'NR') return; end
  if ~isfield(t,'NB') t.NB=size(t.Sab,1); end 
  % NR: Number of interpolation ranges for a particular parameter
  % NB: Total number of beta values given.
  t.head = [ 0,0,0,0, t.NR, t.NB];
  if ~isfield(t,'beta_INT') t.beta_INT=2; end % default: linear/linear
  t.INT  = t.beta_INT;

  t.beta_INT = beta.INT;      % Interpolation schemes
  t.beta     = zeros(1,t.NB);
  t.T        = [];
  t.LT       = t.beta;
  t.NP       = t.NB;
  t.alpha    = [];
  t.Sab      = [];  % S(a,BETAn,Tn)
  % [MAT,7,4/ 0,0,0,0,NR,NB/BetaInt ] TAB2 (interpolation scheme)
  [beta, MF7.lines] = write_endf_TAB2(MF7.lines);
  
  % read alpha,beta,T,Sab values
  for ibeta = 1:t.NB
    % [MAT,7,4/ T0,Beta0,LT,0,NR,NP/AlphaInt/S(a,BETAn,Tn) ] TAB1 Sab(T0)
    [alpha, MF7.lines] = write_endf_TAB1(MF7.lines);
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
    t.Sab(ibeta, :, 1) = alpha.Y; 
    % ENDF: [MAT, 7, 4/ Tn,BETAn,LT,0,NP,0/S(a,BETAn,Tn) ] LIST
    for index=2:(t.LT+1)
      [Sab, MF7.lines] = write_endf_LIST(MF7.lines);
      if isempty(Sab), 
        disp([ mfilename ': ERROR: ' t.description ' LIST/Sab aborted (MF7/MT4)' ]);
        t=[]; return; 
      end % bad format
      if any(Sab.head([ 4 6 ])), 
        disp([ mfilename ': WARNING: ' t.description ' wrong head values in LIST/Sab (MF7/MT4)' ]);
        disp(Sab.head)
      end % bad values LIST
      t.T(index)  = Sab.head(1);% Tn (K)
      t.Sab(ibeta, :, index) = Sab.B;  % append S(E,Tn) as rows
    end
  end % beta loop (NB)
  
  % handle ln(S) storage and compute real Sab
  if t.LLN, t.Sab = exp(t.Sab); end
  
  % read Teff values
  % [MAT,7,4/ 0,0,0,0,NR,NT/Tint/Teff(T) ] TAB1
  Teff_index=6;
  if numel(t.B) >= 13 && t.B(13) == 0, Teff_index=[ Teff_index 13 ]; end
  if numel(t.B) >= 17 && t.B(17) == 0, Teff_index=[ Teff_index 17 ]; end
  t.NT = []; t.Teff_INT = []; t.Teff_T = [];
  for index=1:numel(Teff_index)
    if ~isempty(MF7.lines) && numel(t.B) >= index && t.B(index)
      [Teff, MF7.lines] = write_endf_TAB1(MF7.lines);
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

  % end write_endf_mf7_mt4


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
  

% ------------------------------------------------------------------------------
% ENDF Data Format labels
function MF_description = write_endf_MF(index)
  % get File dscription (MF label)
  persistent MF
  
  if isempty(MF)
  
    MF = cell(1,99); % allow 99 MF sections;
    MF{1}='General information';
    MF{7}='Thermal neutron scattering law data';
  end
  
  if index > 0 && index < numel(MF)
    MF_description = MF{index};
  else
    MF_description = '';
  end
  if isempty(MF_description), MF_description = ''; end

function NLIB_description = write_endf_NLIB(index)
  % get NLIB Library Definition (NLIB label)
  persistent NLIB
  
  if isempty(NLIB)
  
    NLIB = cell(1,50); % allow 50 MF sections;
    NLIB{1+0  } = 'ENDF/B - United States Evaluated Nuclear Data File';
 
  end
  
  if index >= 0 && index < numel(NLIB)
    NLIB_description = NLIB{1+index};
  else
    NLIB_description = '';
  end
  if isempty(NLIB), NLIB_description = ''; end
  
function d = write_endf_NSUB(index)
  % get NSUB sub-library Definition (NSUB label)

  case 12    ; d='Thermal Neutron Scattering Data';

