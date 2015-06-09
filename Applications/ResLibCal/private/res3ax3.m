function [R0, RMS, method] = res3ax3(h,k,l,w,EXP)
% JO 2003 Stoica/Popovici in the Cooper & Nathan approximation

method = 'Stoica/Popovici method in the C. & N. approximation by JO';

% In addition the parameters f, energy pre-multiplier-f*w
% where f=0.48 for meV to ang-2 - and q0 which is the wavevector
% transfer in ang-1, are passed over. 
f  = 0.4826;
CONVERT2=2.072;
pit= pi/60/180; % This is a conversion from minutes of arc to radians.

%----- INPUT SPECTROMETER PARAMETERS.

dm   = EXP.mono.d;            % monochromator d-spacing in Angs.
da   = EXP.ana.d;             % analyser d-spacing in Angs.
etam = EXP.mono.mosaic*pit;   % monochromator mosaic (converted from mins->rads)
etamp= EXP.mono.vmosaic*pit;
etaa = EXP.ana.mosaic*pit;    % analyser mosaic.
etaap= EXP.ana.vmosaic*pit;
etas = EXP.sample.mosaic*pit; % sample mosaic.
etasp= EXP.sample.vmosaic*pit;
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
%-------q0-definition: q=[hkl] w
a = EXP.sample.a;
b = EXP.sample.b;
c = EXP.sample.c;

% a,b,c,alpha,beta,gamma, QH,QK,QL (from ResCal5/rc_re2rc)
[q0,Qmag]= rc_re2rc( [ a b c ], ...
  [ EXP.sample.alpha EXP.sample.beta EXP.sample.gamma ] , ...
  [ h k l ] );  
q0x = q0(1); q0y = q0(2); q0z = q0(3);
q0   = Qmag; % sqrt(q0x^2+q0y^2+q0z^2);
%------Calculate-ki-and-kf,-thetam-and-thetaa
ki=sqrt(kfix^2+(fx-1)*f*w); 
kf=sqrt(kfix^2-(2-fx)*f*w);
%------Test-if-scattering-triangle-is-closed
cos_2theta=(ki^2+kf^2-q0^2)/(2*ki*kf);
if abs(cos_2theta) > 1, 
  disp([ datestr(now) ': ' mfilename ': KI,KF,Q triangle will not close (kinematic equations). Change the value of KFIX,FX,QH,QK or QL.' ]);
  disp([ h k l w ]);
  R0=0; RMS=[];
  return
end
%----------------------------------------
taum=2*pi/da; taua=2*pi/da; q=q0;
thetam=asin(taum/(2*ki))*sm;
thetaa=asin(taua/(2*kf))*sa; 
s2theta=acos( (ki^2+kf^2-q^2)/(2*ki*kf))*ss; %2theta sample
thetas=s2theta/2;
phi=atan2((-kf*sin(s2theta)), (ki-kf*cos(s2theta)));

%--------DEF-DES-MATRICES:
F = diag(1./[etam etamp etaa etaap].^2);

%------matrice-C:
C = zeros(4,8);
C(1,1) = 0.5;
C(1,2) = 0.5;
C(3,5) = 0.5;
C(3,6) = 0.5;
C(2,3) = 1/(2*sin(thetam));
C(2,4) = - C(2,3); % was C(3,3) = - C(2,3); wrong in Popovici paper
C(4,7) = 1/(2*sin(thetaa));
C(4,8) = - C(4,7);

%------matrice-G:
G =diag(1./[alf0 alf1 bet0 bet1 alf2 alf3 bet2 bet3].^2);
    
%------matrice-A:
A = zeros(6,8);
A(1,1) = ki/(2*tan(thetam)); 
A(1,2) = - A(1,1);
A(2,2) = ki; 
A(3,4) = ki; 
A(4,5) = kf/(2*tan(thetaa)); 
A(4,6) = - A(4,5);
A(5,5) = kf; 
A(6,7) = kf; 

%------matrice-B:
B=zeros(4,6);
B(1,1)=cos(phi);
B(1,2)=sin(phi);
B(1,4)=-cos(phi-s2theta);
B(1,5)=-sin(phi-s2theta);
B(2,1)=-B(1,2);
B(2,2)=B(1,1);
B(2,4)=-B(1,5);
B(2,5)=B(1,4);
B(3,3)=1;
B(3,6)=-1;
B(4,1)=2*ki/f;
B(4,4)=-2*kf/f;

% Cooper-Nathans matrices
H    = C'*F*C+G;
Hinv = inv(H);
Ninv = A*Hinv*A';
N    = inv(Ninv);
Minv = B*Ninv*B';
% correction for sample mosaic
Minv(2,2)=Minv(2,2)+q^2*etas^2;
Minv(3,3)=Minv(3,3)+q^2*etasp^2;
M    = inv(Minv);

RMS=8*log(2)*M;

% Calculation of prefactor, normalized to source (Cooper-Nathans)
Rm=ki^3/tan(thetam);
Ra=kf^3/tan(thetaa);
P0=Rm*Ra*(2*pi)^4;           % the det(G) term gets simplified          % Popovici Eq 5
R0=P0/(64*pi^2*sin(thetam)*sin(thetaa));  % Popovici Eq 9

R0=R0*sqrt( det(F)/det(H) ); % Popovici Eq 9

% Transform prefactor to Chesser-Axe normalization
R0=R0/(2*pi)^2*sqrt(det(RMS));
% Include kf/ki part of cross section
R0=R0*kf/ki;
% sample mosaic S. A. Werner & R. Pynn, J. Appl. Phys. 42, 4736, (1971), eq 19
R0=R0/sqrt((1+(q*etasp)^2*RMS(3,3))*(1+(q*etas)^2*RMS(2,2)));

