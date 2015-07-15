function [RMS, S, U]=ResLibCal_RM2RMS(H,K,L,W,EXP,RM)
% RMS=ResLibCal_RM2RMS(H,K,L,W,EXP,RM) rotate Resolution matrix in [Q1 Q2] frame
% [RMS,S]=ResLibCal_RM2RMS(H,K,L,W,EXP,RM)
%   also returns the tranformation matrix in the lattice coordinates
%   Q in frame [ABC] = Q*HKL, e.g. to pass the cloud in the ABC frame

RMS=[]; S=[];
if isempty(RM), return; end

% method: ResCal5 which is the same as legacy RESCAL.
% Transformation matrix 'S' checked against templateTAS.instr and TAS MAD ILL.

HKL = [H K L ]';
%----- Calculate Q2c matrix
[Q2c,Qmag]= rc_re2rc( [ EXP.sample.a EXP.sample.b EXP.sample.c ], ...
  [ EXP.sample.alpha EXP.sample.beta EXP.sample.gamma ] , ...
  HKL );

%----- Now work out transformations
A1=EXP.orient1(:);
A2=EXP.orient2(:);

V1=Q2c*A1;
V2=Q2c*A2;

%----- Form unit vectors V1, V2, V3 in scattering plane

V3=cross(V1,V2);
V2=cross(V3,V1);
V3=V3/norm(V3);
V2=V2/norm(V2);
V1=V1/norm(V1);

U=[V1 V2 V3]';

%----- S transformation matrix from (h,k,l) to V1,V2,V3
S=U*Q2c;  % This is used to bring the Resolution matrix into a defined frame.

%----- Q vector in cartesian coordinates
% Qcart=Q2c*HKL;

%----- Work out angle of Q wrt to V1, V2
TT=S*HKL;
cos_theta=TT(1)/norm(TT);
sin_theta=TT(2)/norm(TT);
R=[cos_theta sin_theta 0; -sin_theta cos_theta 0; 0 0 1];

%----- get and transform resolution matrix in Q frame. RMS is resolution matrix in Qx, Qy & Qz
T=zeros(4,4);
T(4,4)=1;
T(1:3,1:3)=R*S;

%----- Resolution ellipsoid in terms of H,K,L,EN ([Rlu] & [meV]) 
RMS=T'*RM*T;  % 'M' in lattice frame = RMS...


% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% method2: ResLib
% the code extracted from ResLib/ResMatS does not seem to properly compute the 
% RMS matrix, as it assumes scattering to take place in plane, and restricts
% the transformation matrix to [1:2] sub matrix.
% Beware: the code assumes also rows-columns 3 and 4 have been swapped in RM.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


