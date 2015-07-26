function sqw_af_1d
% This particular function calculates the cross section for a gapped 
% excitations in a 1-dimensional antiferromagnet. "a" is the chain axis.
% Arguments H-W and the fields of 'sample' and 'rsample' are vectors, 
% so dont forget to use ".*" instead of "*", etc. This function is
% meant to be used together with the prefactor function PrefDemo.m
%  ResLib v.2.1

% Extract the three parameters contained in "p":
Deltax=p(1);					% Gap at the AF zone-center in meV for x-axis mode
Deltay=p(2);					% Gap at the AF zone-center in meV for y-axis mode
Deltaz=p(3);					% Gap at the AF zone-center in meV for z-axis mode
cc=p(4);						% Bandwidth in meV 
Gamma=p(5);					% Intrinsic HWHM of exccitation in meV
I=p(6);							% Intensity prefactor

% Calculate excitation energy at the (H,K,L) position:
omegax=sqrt(cc^2*(sin(2*pi*H)).^2+Deltax^2);
omegay=sqrt(cc^2*(sin(2*pi*H)).^2+Deltay^2);
omegaz=sqrt(cc^2*(sin(2*pi*H)).^2+Deltaz^2);

% Assume Lorenzian excitation broadening:
lorx=1/pi*Gamma./( (W-omegax).^2+Gamma^2 );
lory=1/pi*Gamma./( (W-omegay).^2+Gamma^2 );
lorz=1/pi*Gamma./( (W-omegaz).^2+Gamma^2 );


% Intensity scales as (1-cos(2*pi*H))/omega0:
sqwx=lorx.*(1-cos(pi*H))./omegax/2;
sqwy=lory.*(1-cos(pi*H))./omegay/2;
sqwz=lorz.*(1-cos(pi*H))./omegaz/2;

% Calculate the polarization factors
% Angle between (0,1,0) [direct space] and [h k l]:
alphay=acos(2*pi*K);
% Angle between (0,0,1) [direct space] and [h k l]:
alphaz=acos(2*pi*L);
% Angle between (1,0,0) [direct space] and [h k l]:
alphax=acos(2*pi*H);

polx=sin(alphax).^2;
poly=sin(alphay).^2;
polz=sin(alphaz).^2;

sqw=sqwy.*poly+sqwz.*polz+sqwx.*polx; %Note: Intensity prefactor applied in PrefDemo.m
