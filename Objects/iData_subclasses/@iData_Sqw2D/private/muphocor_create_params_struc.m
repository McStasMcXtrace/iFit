function [ struc ] = muphocor_create_params_struc( varargin )
%CREATE_PARAMS_STRUC Create the parameter structure for write_mupho
%   
% 
% Create the structure to send the parameters to write_mupho.
% 
% All parameters by default, changed parameters are given in 
% pairs filed, value.
%
% 
% Examples:
% 
%   par = create_params_struc;  % Without parm modification
%   
% % Replacing some values in the structure
%   par = create_params_struc('hox', 60, 'unt', 0.25);
% 
% 
% Parameters:
% -----------
%  
%  input_mupho.txt example and parameters:
%  
% !     E0:INC. ENERGY(MEV)             FP:FLIGHTPATH(CM)
% !     XNEL:CHANNEL OF ELASTIC LINE    CW:CHANNEL WIDTH(MYS)
% !     FIMI:MIN.SCATT.ANGLE            FIMA:MAX.SCATT.ANGLE
% !     UNT:CONST.BACKGROUND            ABK: Absorption coefficient
% !     AEMP: Coeff. for calculating the detector efficiency
% 
%  e0,     fp,         cw,      xnel,     fimi,    fima,       abk,      aemp
%  4.77    248.30      3.91    794.08     10.33    115.05      0.00      3.30
% 
% !       Correction OF constant backgroud FOR counter efficiency
%   bac=unt/(1.0-EXP(-aemp/SQRT(e0-hot(n))))
% 
% !     NSPEC: NUMBER OF DIFFERENT ATOMIC SPECIES
% !     IF IQUO is unequal 0 an external model FOR the ratio OF F/G has
% !     to be provided. The model is given through the array PARQ
%  lm,     nspec,  iquo,lpar
%  500     1       0    1
% 
% !     HOX:   UPPER LIMIT OF DENSITY OF STATES
% !     TEMPO: TEMPERATURE IN KELVIN
% !     DW:    MEAN DEBY WALLER COEFFICIENT
% 
% hox,    temp0,      dw
% 40.00    298.52      0.03
% 
% 
% !     AMASI: ATOMIC MASS
% !     SIGI : SIGMA
% !     CONCI: CONCENTRATION
% !     ALFI:  SCATTERING POWER
% DO  i=1,nspec <-- 1 in our CASE so far
% READ 40,amasi(i),sigi(i),conci(i)
% END DO
% amasi(i),  sigi(i),  conci(i)
% 44.80      3.90      1.00
% 
% !     NPHO:     NUMBER OF MULTI PHONON TERMS
% !     ITM:      TOTAL NUMBER OF ITERATIONS
% !     IVIT=0(1):ITERATION BY DIFFERENCE(QUOTIENT) METHOD
% !     IDW=0:    DW COEFF. KEPT CONST.=INPUT VALUE,DW=1:DW COEFF.ITERATED
% !     IEMP=1(0):CORRECTIONS FOR COUNTER EFFICIENCY WILL BE(NOT BE) DONE
% !     IRES=1(0):DATA WILL BE(NOT BE) CORRECTED FOR SPECTR. RESOLUTION
% !     IPR=1(0): Convolutions integrals o the multiphonon term are printed
% !     ILOSS=1:  Analysis OF the imcomplete energy loss spectrum will be done
% 
% npho,itm, ivit, idw, iemp, ires, ipr,iloss
% 5    10   10    1    1     0     0   0
% 
% !     NU:  First channel OF spektrum NO: Last channel OF spektrum
% !     NUU: First channel OF TOF distribution used FOR calculation
% !     NOO: Last channel  OF TOF distribution used FOR calculation
% !     IGLU: Number OF smoothing processes in CASE OF a time-dependant background
% 
% isour,nuu, noo, nu,  no,    iglu
% 0     200  750  200  750    0
% 
% 
% !     UNT: CONSTANT BACKGROUND
% !     FUN: MULTIPLICATION FACTOR FOR TIME DEPENDENT BACKGROUND
% 
% unt,      fun
% 1.0000    0.0000
% 
% 
%      Z00: TIME OF FLIGHT DISTRIBUTION
%      UN0: TIME DEPENDENT BACKGROUND (if fun = 1)
% 
% 

% 
% 
%SEE ALSO WRITE_MUPHO
%narginchk(0,100);


% Default structure
% ----------------------

% First line
struc.e0 = 4.77;                    % E0:INC. ENERGY(MEV) 
struc.fp = 248.30;                  % FP:FLIGHTPATH(CM) by default here: in6
struc.cw = 3.91;                    % CW:CHANNEL WIDTH(MYS)
struc.xnel = 794.08;                % XNEL:CHANNEL OF ELASTIC LINE 
struc.fimi = 0.0;                   % FIMI:MIN.SCATT.ANGLE (to be filled with data or change params)
struc.fima = 180.0;                 % FIMA:MAX.SCATT.ANGLE
struc.abk  = 0.0;                   % ABK: Absorption coefficient
struc.aemp = 3.30;                  % AEMP: Coeff. for calculating the detector efficiency

% Second line
struc.lm = 500;                     % LM: NUMBER OF DISCRET ENERGY INTERVALS         
struc.nspec = 1;                    % NSPEC: NUMBER OF DIFFERENT ATOMIC SPECIES (?)
struc.iquo  = 0;                    % If IQUO is unequal 0 an external model for the ratio of F/G has to be provided.
struc.lpar  = 1;


% Third line
struc.hox  = 40.0;                  % HOX:   UPPER LIMIT OF DENSITY OF STATES
struc.temp0= 295.00;                % TEMPO: TEMPERATURE IN KELVIN (a default value provided to avoid crashes)
struc.dw   = 0.0;                   % DW:    MEAN DEBY WALLER COEFFICIENT

% Fourth line
struc.amasi = 44.80;                % AMASI: ATOMIC MASS
struc.sigi  = 3.90;                 % SIGI : SIGMA
struc.conci = 1.00;                 % CONC:  CONCENTRATION

% Fifht line
struc.npho  = 5;                    % NPHO: NUMBER OF MULTI PHONON TERMS,ITM:NUMBER OF ITERATIONS            
struc.itm   = 10;                   % ITM: TOTAL NUMBER OF ITERATIONS
struc.ivit  = 10;                   % IVIT=0(1):ITERATION BY DIFFERENCE(QUOTIENT) METHOD
struc.idw   = 1;                    % IDW=0:    DW COEFF. KEPT CONST.=INPUT VALUE,DW=1:DW COEFF.ITERATED
struc.iemp  = 0;                    % IEMP=1(0):CORRECTIONS FOR COUNTER EFFICIENCY WILL BE(NOT BE) DONE
struc.ires  = 0;                    % IRES=1(0):DATA WILL BE(NOT BE) CORRECTED FOR SPECTR. RESOLUTION
struc.ipr   = 0;                    % IPR=1(0): CONVOLUTION INTEGRAL OF THE MULTIPHONON TERM PRINTED OUT
struc.iloss = 0;                    % ILOSS=1(0): ANALYSIS OF OF INCOMPLETE ENERGY LOSS SPECTRUM WILL BE DONE

% Sixth line
struc.isour = 0;                    % NU:  First channel of spektrum 
struc.nuu   = 200;                  % NUU: First channel of TOF distribution used for calculation     
struc.noo   = 750;                  % NOO: Last channel  of TOF distribution used for calculation
struc.nu    = 200;                  % NU:  First channel of spektrum NO: Last channel of spektrum
struc.nuu   = 200;                  % NUU: First channel of TOF distribution used for calculation
struc.no    = 750;                  % NO: Last channel of spektrum
struc.iglu  = 0;

% Seventh line
struc.unt = 0.0;                    % UNT: CONSTANT BACKGROUND
struc.fun = 0.0;                    % FUN: MULTIPLICATION FACTOR FOR TIME DEPENDENT BACKGROUND

% Substitute values for the given parameters
if nargin~=0   
    labels = {varargin{1:2:end-1}};
    values = {varargin{2:2:end}};
    for k=1:numel(values)
        struc.(labels{k})=values{k};
    end
end



% Next are the sin theta reduced time-of-flight data. Z00 and UN0








end

