function [R0,RMS,vi,vf,Error]=rc_cnmat(f,q0,p,mon_flag)
%
% RESCAL  function to calculate the resolution matrix NP in terms of 
% (DQX,DQY,DQZ,DW) defined along the wavevector transfer Q direction 
% in a right hand coordinate system.
%
% The resolution matrix agrees with that calculated using RESCAL from
% the ILL.
%
% Notes: 1. The sign errors in Mitchell's paper have been corrected.
%        2. Error =1 on exit if the scattering triangle does not
%	    close, or imag(NP)~=0.  
%        3. We have followed Dorner and included the kf/ki factor in the
%           normalisation
%        4. The monitor efficiency 1/ki is also included 
%
% Input :       f = Converts from Angs^1 to energy units
%              q0 = Q vector in Angs^-1
%               p = Spectrometer and scan parameters
%         mon_flag= Monitor flag. If mon_flag=1, experiment performed at const. monitor
%                                 If mon_flag=0,    "          "          "     time
%
% Output:       R0= resoltuion volume calculated using Dorner's method.
%               NP= resolution matrix using a matrix method like Mitchell's,
%               vi= incident resolution volume
%               vf= final resolution volume
%            Error= 1, Scattering triangle will not close
%
% ResCal5/AT and DFM, 29.11.95
%

% specificity: original CN only takes isotropic mosaicity

pit=0.0002908882; % This is a conversion from minutes of arc to radians.

%----- INPUT SPECTROMETER PARAMETERS.

if isstruct(p)
  EXP=p;
  % Cooper-Nathans parameters
  dm   = EXP.mono.d;            % monochromator d-spacing in Angs.
  da   = EXP.ana.d;             % analyser d-spacing in Angs.
  etam = EXP.mono.mosaic*pit;   % monochromator mosaic (converted from mins->rads)
  etaa = EXP.ana.mosaic*pit;    % analyser mosaic.
  etas = EXP.sample.mosaic*pit; % sample mosaic.
  sm   = EXP.mono.dir;          % scattering sense of monochromator (left=+1,right=-1)
  ss   = EXP.sample.dir;        % scattering sense of sample (left=+1,right=-1)
  sa   = EXP.ana.dir;           % scattering sense of analyser (left=+1,right=-1)
  kfix = EXP.Kfixed;            % fixed momentum component in ang-1.
  fx   = 2*(EXP.infin==-1)+(EXP.infin==1);             % fx=1 for fixed incident and 2 for scattered wavevector.
  alf0 = EXP.hcol(1)*pit;       % horizontal pre-monochromator collimation.
  alf1 = EXP.hcol(2)*pit;       % horizontal pre-sample collimation.
  alf2 = EXP.hcol(3)*pit;       % horizontal post-sample collimation.
  alf3 = EXP.hcol(4)*pit;       % horizontal post-analyser collimation.
  bet0 = EXP.vcol(1)*pit;       % vertical pre-monochromator collimation.
  bet1 = EXP.vcol(2)*pit;       % vertical pre-sample collimation.
  bet2 = EXP.vcol(3)*pit;       % vertical post-sample collimation.
  bet3 = EXP.vcol(4)*pit;       % vertical post-analyser collimation.
  w    = EXP.W;
  
  % a,b,c,alpha,beta,gamma, QH,QK,QL (from ResCal5/rc_re2rc)
  [q2c,q0]= rc_re2rc( [ EXP.sample.a EXP.sample.b EXP.sample.c ], ...
    [ EXP.sample.alpha EXP.sample.beta EXP.sample.gamma ] , ...
    [ EXP.QH EXP.QK EXP.QL ] );  
elseif isvector(p)
  dm=p(1);            % monochromator d-spacing in Angs.
  da=p(2);            % analyser d-spacing in Angs.
  etam=p(3)*pit;      % monochromator mosaic (converted from mins->rads)
  etaa=p(4)*pit;      % analyser mosaic.
  etas=p(5)*pit;      % sample mosaic.
  sm=p(6);            % scattering sense of monochromator (left=+1,right=-1)
  ss=p(7);            % scattering sense of sample (left=+1,right=-1)
  sa=p(8);            % scattering sense of analyser (left=+1,right=-1)
  kfix=p(9);          % fixed momentum component in ang-1.
  fx=p(10);           % fx=1 for fixed incident and 2 for scattered wavevector.
  alf0=p(11)*pit;     % horizontal pre-monochromator collimation.
  alf1=p(12)*pit;     % horizontal pre-sample collimation.
  alf2=p(13)*pit;     % horizontal post-sample collimation.
  alf3=p(14)*pit;     % horizontal post-analyser collimation.
  bet0=p(15)*pit;     % vertical pre-monochromator collimation.
  bet1=p(16)*pit;     % vertical pre-sample collimation.
  bet2=p(17)*pit;     % vertical post-sample collimation.
  bet3=p(18)*pit;     % vertical post-analyser collimation.
  w=p(34);            % energy transfer.
end

% In addition the parameters f, energy pre-multiplier-f*w
% where f=0.48 for meV to ang-2 - and q0 which is the wavevector
% transfer in ang-1, are passed over. 

% Calculate ki and kf, thetam and thetaa

ki=sqrt(kfix^2+(fx-1)*f*w);  % kinematical equations.
kf=sqrt(kfix^2-(2-fx)*f*w);

% Test if scattering triangle is closed

cos_2theta=(ki^2+kf^2-q0^2)/(2*ki*kf);
if abs(cos_2theta) <= 1, Error=0; else
  disp([ datestr(now) ': ' mfilename ': Can not close triangle (kinematical equations). ' ])
  if isstruct(p), disp([ EXP.QH EXP.QK EXP.QL EXP.W ]);
  else            disp(p(31:34));
  end
  R0=0; RMS=[];
  return
end

thetaa=asin(pi/(da*kf));      % theta angles for analyser
thetam=asin(pi/(dm*ki));      % and monochromator.
thetas=acos(cos_2theta)*ss/2;

N=zeros(6,6);	% this is matrix N = inv(A*Hinv*A') with Hinv=inv(C'*F*C+G) in Popovici nomenclature
B=N;

% Fill up the horizontal components first.

pm=1/(ki*etam)   *[sm*tan(thetam)    1]; % [ a1 a2 ]
palf0=1/(ki*alf0)*[2*sm*tan(thetam)  1]; % [ a7 a8 ]
palf1=1/(ki*alf1)*[0                 1]; % [ 0  a3 ]

pa=1/(kf*etaa)   *[-sa*tan(thetaa)   1]; % [ -a5 -a6 ]
palf3=1/(kf*alf3)*[-2*sa*tan(thetaa) 1]; % [ -a9  a10 ]
palf2=1/(kf*alf2)*[0                 1]; % [ 0    a4 ]

N(1:2,1:2)=pm'*pm+palf0'*palf0+palf1'*palf1;
N(4:5,4:5)=pa'*pa+palf3'*palf3+palf2'*palf2;

% Now fill up the vertical components.

b1=1/(bet1^2)+1/((2*sin(thetam)*etam)^2+bet0^2);  % these are Dorner's 
b2=1/(bet2^2)+1/((2*sin(thetaa)*etaa)^2+bet3^2);  % corrected formulae.

N(3,3)=1/(ki^2)*b1;
N(6,6)=1/(kf^2)*b2;

% The resolution matrix in terms of the incident and scattered
% momentum 3-vectors has been constructed.

% Now calculate the transformation matrix.
%

ang1=acos(-(kf^2-q0^2-ki^2)/(2*q0*ki));       % angle between ki and q0
ang2=pi-acos(-(ki*ki-q0*q0-kf*kf)/(2*q0*kf)); % angle between kf and q0

TI=[ cos(ang1) -ss*sin(ang1) ; ss*sin(ang1) cos(ang1) ]; % transform kix,kiy
                                                         % to DQX and DQY.
TF=[ cos(ang2) -ss*sin(ang2) ; ss*sin(ang2) cos(ang2) ]; % transform Kfx,kfy
                                                         % to DQX and DQY.
% this is kind of matrix 'B'
B(1:2,1:2)=TI;
B(1:2,4:5)=-TF;
B(3,3)=1;
B(3,6)=-1;
B(4,1)=2*ki/f;
B(4,4)=-2*kf/f;
B(5,1)=1;
B(6,3)=1;

% estimate from old ResCal5 code only (pure ResCal)
Uold = B;
Vold = inv(Uold);
Mold = N; % N = inv(A*Hinv*A')
Nold = Vold'*Mold*Vold;

[dummy,Nold]=rc_int(6,1,Nold);        % integrate over kiz giving a 5x5 matrix
[dummy,Nold]=rc_int(5,1,Nold);        % integrate over kix giving a 4x4 matrix
NP=Nold-Nold(1:4,2)*Nold(1:4,2)'/(1/((etas*q0)^2)+Nold(2,2));
NP(3,3)=Nold(3,3);
NP=8*log(2)*NP;                        % Correction factor as input parameters
                                       % are expressed as FWHM.

RMS = NP;

%----- Normalisation factor


% Calculation of prefactor, normalized to source (Cooper-Nathans)
% we provide below a set of estimation methods

%% the original ResCal5/cn_mat intensity estimate
if 1
  vi=ki^3*cot(thetam);
  vi=vi/sqrt((2*sin(thetam)*etam)^2+bet0^2+bet1^2);
  vf=kf^3*cot(thetaa);
  vf=vf/sqrt((2*sin(thetaa)*etaa)^2+bet2^2+bet3^2);
  R0=vi*vf*(2*pi)^4;

  sqrt_detF_over_detH = (15.75*bet0*bet1*etam*alf0*alf1)/sqrt(alf0^2+alf1^2+4*etam^2) ...
                      * (15.75*bet2*bet3*etaa*alf2*alf3)/sqrt(alf2^2+alf3^2+4*etaa^2);

  R0 = R0 * sqrt_detF_over_detH;
end

%% legacy computation from Chesser and Axe Acta Cryst 1972 (provides same result as AFILL)
if 0
  a11 = b1/ki^2;
  a12 = b2/kf^2;

  a1= pm(1); a2=pm(2);  a3=palf1(2); a4= palf2(2);
  a5=-pa(1); a6=-pa(2); a7=palf0(1); a8= palf0(2);
  a9=-palf3(1); a10=-palf3(2);

  b0=a1*a2+a7*a8;
  b1=a2*a2+a3*a3+a8*a8;
  b2=a4*a4+a6*a6+a10*a10;
  b3=a5*a5+a9*a9;
  b4=a5*a6+a9*a10;
  b5=a1*a1+a7*a7;

  ALAM=ki/kf;
  AL=sin(2*thetas)*ss*sm; BE=cos(2*thetas);

  C=-1.*(ALAM-BE)/AL;
  E=-1.*(BE*ALAM-1.)/AL;
  AP=2.*b0*C+b1*C*C+b2*E*E+b3*ALAM*ALAM+2.*b4*ALAM*E+b5; % A' in CN paper

  R0 = 2*pi/(ki^2*kf^3*AL) ...
		  /sqrt(AP*(a11+a12)) ...
		  *sqrt( ...
		     bet0^2/(bet0^2+(2*etam*sin(thetam))^2) ...
		    *bet3^2/(bet3^2+(2*etaa*sin(thetaa))^2));
  R0 = abs(R0);
end

% Transform prefactor to Chesser-Axe normalization
R0=R0/(2*pi)^2*sqrt(det(RMS));
% Include kf/ki part of cross section
R0=R0*kf/ki;
% sample mosaic S. A. Werner & R. Pynn, J. Appl. Phys. 42, 4736, (1971), eq 19
R0=R0/sqrt((1+(q0*etas)^2*RMS(3,3))*(1+(q0*etas)^2*RMS(2,2)));

%----- Final error check

if imag(RMS) == 0; Error=0; else; Error=1; end

return
