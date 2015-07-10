function RMS=ResLibCal_RM2RMS(H,K,L,W,EXP,RM)
% RMS=ResLibCal_RM2RMS(H,K,L,W,EXP,RM) rotate Resolution matrix in [Q1 Q2] frame

if isempty(RM), RMS=RM; return; end

% method1: ResLib
% code extracted from ResLib/ResMatS
% Calls: CleanArgs, StandardSystem, modvec, scalar

if 0
  % disp('Using ResLib RM2RMS')
  [len,H,K,L,W,EXP]=CleanArgs(H,K,L,W,EXP);
  [x,y,z,sample,rsample]=StandardSystem(EXP);

  Q=modvec(H,K,L,rsample);
  uq(1,:)=H./Q;  % Unit vector along Q
  uq(2,:)=K./Q;
  uq(3,:)=L./Q;

  xq=scalar(x(1,:),x(2,:),x(3,:),uq(1,:),uq(2,:),uq(3,:),rsample);
  yq=scalar(y(1,:),y(2,:),y(3,:),uq(1,:),uq(2,:),uq(3,:),rsample);
  zq=0;  %scattering vector assumed to be in (orient1,orient2) plane;

  tmat=zeros(4,4,len); %Coordinate transformation matrix
  tmat(4,4,:)=1;
  tmat(3,3,:)=1;
  tmat(1,1,:)=xq;
  tmat(1,2,:)=yq;
  tmat(2,2,:)=xq;
  tmat(2,1,:)=-yq;

  RMS_swapped=zeros(4,4,len);
  
  % this code applies on RM with columns/rows 3-4 swapped.
  RM_swapped = RM;
  RM_swapped([3 4], :) = RM([4 3],:);
  RM_swapped(:, [3 4]) = RM(:, [4 3]);
  
  for i=1:len
     RMS_swapped(:,:,i)=(tmat(:,:,i))'*RM_swapped(:,:,i)*tmat(:,:,i);
  end
  
  RMS = RMS_swapped;
  RMS([3 4], :) = RMS_swapped([4 3],:);
  RMS(:, [3 4]) = RMS_swapped(:, [4 3]);
end

% method2: ResCal5 which is the same as legacu RESCAL
if 1
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
  RMS=T'*RM*T;  % M should be RMS...
end
