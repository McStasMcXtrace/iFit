function out = ResLibCal_AxesResMat
% ResLibCal_AxesResMat: creates an axis system of Monte-Carlo points
%   which represent the resolution function
%
%   EXP: structure containing application configuration. When not given, extracts
%        it from the main window.
%
% Returns:
%   out: full output from ResLibCal, with resolution 
%        and 4D axes as { H,K,L,W } Monte-Carlo cloud in out.resolution

% insprired from ResCal5/mc_conv

ax=[];
out        = ResLibCal_Compute;
EXP        = out.EXP;
resolution = out.resolution;

% handle case with vector of HKLE locations
if iscell(out.resolution)
  for index=1:numel(out.resolution)
    out.resolution{index}.cloud = ResLibCal_AxesResMatSingle(out, EXP, out.resolution{index});
  end
else
  out.resolution.cloud = ResLibCal_AxesResMatSingle(out, EXP, resolution);
end

% ==============================================================================
function ax = ResLibCal_AxesResMatSingle(out, EXP, resolution)
% compute a Monte-Carlo cloud for a single HKLE location, using the resolution function

% resolution.HKLE is the location of the scan step where the convolution must be evaluated

%----- Calculate Q2c matrix
[Q2c,Qmag]= rc_re2rc( [ EXP.sample.a EXP.sample.b EXP.sample.c ], ...
  [ EXP.sample.alpha EXP.sample.beta EXP.sample.gamma ] , ...
  resolution.HKLE(1:3)' );

%----- Now work out transformations
A1=EXP.orient1(:);
A2=EXP.orient2(:);

V1=Q2c*A1;
V2=Q2c*A2;

%----- Form unit vectors V1, V2, V3 in scattering plane

V3=cross(V1,V2);
V2=cross(V3,V1);
V3=V3/sqrt(sum(V3.*V3));
V2=V2/sqrt(sum(V2.*V2));
V1=V1/sqrt(sum(V1.*V1));

U=[V1';V2';V3'];

%----- S transformation matrix from (h,k,l) to V1,V2,V3
S=U*Q2c;  % This is used to bring the CN matrix into a defined frame.

%----- Q vector in cartesian coordinates
Qcart=Q2c*resolution.HKLE(1:3)';
Qmag=sqrt(sum(Qcart.*Qcart));

%----- Work out angle of Q wrt to V1, V2
TT=S*resolution.HKLE(1:3)';
cos_theta=TT(1)/sqrt(sum(TT.*TT));
sin_theta=TT(2)/sqrt(sum(TT.*TT));
R=[cos_theta sin_theta 0; -sin_theta cos_theta 0; 0 0 1];

%----- get and transform resolution matrix in Q frame. RMS is resolution matrix in Qx, Qy & Qz
T=zeros(4,4);
T(4,4)=1;
T(1:3,1:3)=R*S;

%----- Resolution ellipsoid in terms of H,K,L,EN ([Rlu] & [meV]) 
M=T'*resolution.RM*T; 
[V,E]=eig(M);
sigma=1./sqrt(diag(E)); % length along principal axes of gaussian

% compute MC points on axes (with NMC points)
NMC  = 10000;
xp   = bsxfun(@times,sigma,rand(4,NMC));
XMC  = inv(V)'*xp; % this is delta(HKLW) as a gaussian distribution
clear xp
HKLW = bsxfun(@plus,resolution.HKLE',XMC); % add HKLE location to delta

ax = [ mat2cell(HKLW,[1 1 1 1]) ]; % get 1D arrays per axis

% and now we can evaluate the function onto the axes.... and sum all values
% sum(feval(model, parameters, ax{:}))*resolution.R0/NMC
