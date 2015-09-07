function res = ResLibCal_RM2RMS(EXP, res)
% RMS=ResLibCal_RM2RMS(H,K,L,W,EXP,RM) rotate Resolution matrix in [Q1 Q2] frame
% [RMS,S]=ResLibCal_RM2RMS(H,K,L,W,EXP,RM)
%
% Compute for each coordinate frame:
%   cart2frame: conversion matrix from cartesian lattice [Angs-1] to current frame in [uni...
%   rlu2frame:  conversion matrix from HKL lattice [rlu] to current frame in [unit] below
%   Q:          location of the HKL evaluation in the frame, in [unit] below
%   RM:         resolution matrix in the frame, in [unit]^3.meV
%   README:     a short description of the coordinate frame
%   unit:       the unit used for wavevector representation
%   frame:      the unit vectors defining the frame, expressed in the [frameUnit] below
%   frameUnit:  the wavevector unit used for defining the [frame]
%   frameStr:   the human readable expression of the [frame] unit vectors, in [frameUnit]
%   Bragg:      the resolution FWHM along [x,y,z,E,Vana] in [unit]^3.meV
% The available coordinate frames are:
%      rlu: [R] rlu lattice frame, xyz=[a* b* c*]
%     cart: [B] cartesian lattice frame, ortho-normal, x=a*/|a*|
%     spec: [Q] spectrometer frame, ortho-normal, x=Q/|Q|
%  rlu_ABC: [V] rlu rotated lattice frame xyz=[A B C]
%      ABC: [U] cartesian rotated lattice frame, ortho-normal, x=A/|A|
%
% The coordinate frames used by the legacy RESCAL are [R] (rlu) and [Q] (spec).
% The [U] (ABC) coordinate frame is needed to compute the TAS angles.

% require the resolution matrix in the spectrometer frame to be available
% this one is computed in ResLibCal_ComputeResMat (defines res.xyz.RM).
if isempty(res.spec.RM), return; end
RM_Q = res.spec.RM;

% compute the resolution information for all desired frames

% lattice reciprocal frame in [rlu]
% the B matrix is the reciprocal lattice vector frame, expressed in cartesian.
%   W.R. Busing and H.A. Levy, Acta Cryst. 22 , 457 (1967)

HKL = res.HKLE(1:3)';

% ==============================================================================
% Theory:
% transformation matrix P: the columns of the transformation matrix =
%    coordinates of the new frame vectors, expressed in the old frame.
% then
%    [old vector] = P*[new vector]
% so the nomenclature is opposite to the usual way. i.e. P looks like it 
%   converts from [new] to [old].
% e.g. if [new vector] is [1 0 0], we get the first basis vector of the new 
%  basis, expressed in the old basis. P: new -> old
% ==============================================================================

% [B] = [a* b* c*] is the lattice frame, expressed in cartesian coords
%     This frame defines [rlu] units.
%     aka "cartesian lattice x//(a*/|a*|)" frame
% [B] == rlu2cartesian = R2B
%
% This frame provides the "cartesian" [HKL] measurement location in the crystal
%   lattice frame. This is usefull if e.g. the fit model is in [Angs-1] with
%   principal axes matching those of the crystal.
[B, QM, Qcart]= rc_re2rc( [ EXP.sample.a EXP.sample.b EXP.sample.c ], ...
  [ EXP.sample.alpha EXP.sample.beta EXP.sample.gamma ] , ...
  HKL );
  
% [R] is the lattice frame, expressed in rlu coords, i.e. identity
%     aka "rlu lattice xyz=[a* b* c*]" frame
% [R] == rlu2rlu = R2R
%
% This frame provides the "rlu" [HKL] measurement location in the crystal
%   lattice frame. This is usefull if e.g. the fit model is in [rlu] with
%   principal axes matching those of the crystal.
R = eye(3);
  
% [U] = [A B C] is a lattice frame, oriented along user defined vectors [ABC]
%     expressed in cartesian coords. This frame is built ortho-normal.
%     aka "cartesian rotated lattice x//(A/|A|)" frame
% [U] == ABC2cartesian == U2B
%
% This frame provides the "cartesian" [HKL] measurement location in the crystal
%   lattice frame, rotated to match the user defined A and B scattering vectors. 
% This is usefull if e.g. the fit model is in [Angs-1] and must be rotated.
% When A=a*, this roughtly corresponds with the [B] frame, except this one is 
%   orthonormal, and [B] is usually not (except cubic...).
U1=B*EXP.orient1(:);  % user defined 'A' vector converted from rlu to Angs-1, in [a* b* c*] frame
U2=B*EXP.orient2(:);  % user defined 'B' vector converted from rlu to Angs-1, in [a* b* c*] frame

U3=cross(U1,U2);  % V3 is perp(A,B), normalised
U2=cross(U3,U1);
U3=U3/norm(U3);
U2=U2/norm(U2);
U1=U1/norm(U1); % V1 is along user defined vector 'A', normalised (cartesian)

U=[U1(:) U2(:) U3(:)];  % ortho-normal ABC vector basis in [a* b* c*] frame, in Angs-1
                % V1//A, V3=(AxB)
                
% [Q] is the cartesian "spectrometer" frame. This frame is built ortho-normal.
%     aka "spectrometer x//Q" frame
% [Q] == Qyz2cartesian = Q2B
%
% This frame provides the "cartesian" [HKL] measurement location in the 
%   spectrometer frame. 
% This is usefull for e.g. isotropic materials (powders, gas, liquids, glasses...)
x = Qcart/norm(Qcart);  % normalised frame. x // Q
z = U3; % a 'vertical' axis
y = cross(z,x); y=y/norm(y);
z = cross(x,y); z=z/norm(z);
Q = [ x(:) y(:) z(:) ];

% [V] is a rlu lattice frame, oriented along user defined vectors [ABC]
%     aka "rlu rotated lattice xyz=[A B C]" frame. NOT ortho-normal.
% [V] == rluABC2rlu = V2R
%
% This frame provides the "rlu" [HKL] measurement location in the crystal
%   lattice frame, rotated to match the user defined A and B scattering vectors. 
V1 = EXP.orient1(:);  % user defined 'A' vector, as is, in [rlu]
V2 = EXP.orient2(:);  % user defined 'B' vector, as is, in [rlu]
V3 = inv(B)*U3;       % use the [U]=[ABC] cartesian vertical vector, and convert it back into [rlu]

V = [ V1(:) V2(:) V3(:)];

% compute scattering vector in all these bases: 1st the RLU frames, then the cartesian [Angs-1]
Q_R = HKL;          % == R*HKL; in units of [a* b* c*]  R=rlu2rlu
Q_V = inv(V)*HKL;   %           in units of [ABC]       V=rluABC2rlu,     inv(V)=rlu2rluABC

Q_B = Qcart;        % == B*HKL;                         B=rlu2cart
Q_U = U'*Qcart;     %                                   U=ABC2cart,       inv(U)=U'=cart2ABC
Q_Q = Q'*Qcart;     %                                   Q=Qyz2cartesian,  inv(Q)=Q'=cart2Qyz

% transformation matrices from HKL to given frames, and inverse.
U2B = U;  % ortho-normal
V2R = V;
Q2B = Q;  % ortho-normal

% 'hkl':'rlu' == [R], 'cart':'ABC' == [B], 'spec':'xyz' == [Q]
R2V = inv(V); % hkl to 'rlu [A B C]'          inv(V2R)
R2B = B;      % hkl to "cartesian [a* b* c*]" 
R2U = U'*B;   % hkl to 'cartesian [A B C ]'   R->B->U: B2U*R2B=inv(U)*B=U'*B
R2Q = Q'*B;   % hkl to "spectrometer Q"       R->B->Q: B2Q*R2B=inv(Q)*B=Q'*B

% compute the resolution matrix in all these bases
% the resolution mtrix is given in the [Q] frame, in cartesian [Angs-1].
% to convert it to an other [frame]:
%      T = [frame]2Q
%    then RM[frame] = T'*RM_Q*T

% from Q[xyz] to R[rlu]:    T=R2Q=B2Q*R2B=Q'*B
T=diag([0 0 0 1]); T(1:3,1:3)=R2Q;      RM_R = T'*RM_Q*T;
% from Q[xyz] to B[cart]:   T=B2Q=Q'
T=diag([0 0 0 1]); T(1:3,1:3)=Q';       RM_B = T'*RM_Q*T;
% from Q[xyz] to V[rluABC]: T=V2Q=R2Q*V2R=(Q'*B)*V
T=diag([0 0 0 1]); T(1:3,1:3)=R2Q*V2R;  RM_V = T'*RM_Q*T;
% from Q[xyz] to U[ABC]:    T=U2Q=B2Q*U2B
T=diag([0 0 0 1]); T(1:3,1:3)=Q'*U;     RM_U = T'*RM_Q*T;

% assemble results into a comprehensive structure
Angs = [ '1/' char(197) ];
% [R] basis
res.rlu.cart2frame    = inv(B); % B2R
res.rlu.rlu2frame     = eye(3); % R2R
res.rlu.Q             = Q_R;    % HKL
res.rlu.RM            = RM_R;   % << THIS is the [RLU] 'RMS' matrix computed by RESCAL in rlu^3.meV
res.rlu.README        = '[R] rlu lattice frame, xyz=[a* b* c*]';
res.rlu.unit          = 'rlu';
res.rlu.frame         = B;      % these vectors are in Angs-1, in the lattice frame
res.rlu.frameUnit     = Angs; % 'Angs-1';
res.rlu.frameStr      = {'a*','b*','c*'};
res.rlu.Bragg         = rc_bragg(res.rlu.RM);

% [B] basis
res.cart.cart2frame   = eye(3); % B2B
res.cart.rlu2frame    = B;      % R2B
res.cart.Q            = Q_B;    % Qcart = B*HKL
res.cart.RM           = RM_B;   % this matrix is the resolution in the lattice frame, expressed in Angs-1
res.cart.README       = '[B] cartesian lattice frame, ortho-normal, x=a*/|a*|';
res.cart.unit         = [ '1/' char(197) ];
res.cart.frame        = eye(3); % in Angs-1 along axes
res.cart.frameUnit    = Angs;
res.cart.frameStr     = {'a*/|a*|', [],[]};
res.cart.Bragg        = rc_bragg(res.cart.RM);

% [Q] basis
res.spec.cart2frame   = Q';     % B2Q
res.spec.rlu2frame    = R2Q;    % R2Q
res.spec.Q            = Q_Q;    % this is Q in the spectrometer space [longitudinal, transverse, vertical]
res.spec.RM           = RM_Q;   % << THIS is the spectrometer 'RM' matrix computed by RESCAL in Angs-3.meV
res.spec.README       = '[Q] spectrometer frame, ortho-normal, x=Q/|Q|';
res.spec.unit         = Angs;
res.spec.frame        = Q;      % in Angs-1 along axes, in the lattice frame
res.spec.frameUnit    = Angs;
res.spec.frameStr     = {'Q//', 'Q⊥', 'Q↑'};
res.spec.Bragg        = rc_bragg(res.spec.RM);

% [V] basis
res.rlu_ABC.cart2frame= R2V*inv(B);    % B2V=R2V*B2R
res.rlu_ABC.rlu2frame = R2V;           % R2V
res.rlu_ABC.Q         = Q_V;
res.rlu_ABC.RM        = RM_V;   % this is the [RLU] 'RMS' matrix rotated in the ABC frame
res.rlu_ABC.README    = '[V] rlu rotated lattice frame xyz=[A B C]';
res.rlu_ABC.unit      = 'rlu';
res.rlu_ABC.frame     = V;      % in rlu along axes
res.rlu_ABC.frameUnit = 'rlu';
res.rlu_ABC.frameStr  = {'A','B','C'};
res.rlu_ABC.Bragg     = rc_bragg(res.rlu_ABC.RM);

% [U] basis
res.ABC.cart2frame    = U';     % B2U
res.ABC.rlu2frame     = R2U;    % R2U
res.ABC.Q             = Q_U;
res.ABC.RM            = RM_U;   % this is the [Angs-1] 'RMS' matrix rotated in the ABC frame
res.ABC.README        = '[U] cartesian rotated lattice frame, ortho-normal, x=A/|A|';
res.ABC.unit          = Angs;
res.ABC.frame         = U;      % in Angs-1 along axes
res.ABC.frameUnit     = Angs;
res.ABC.frameStr      = {'A/|A|',[],[]};
res.ABC.Bragg         = rc_bragg(res.ABC.RM);

res.QM                = QM;
res.README            = {'This structure holds the resolution function expressed in different frames', ...
  'Each structure contains the following fields:', ...
  '  cart2frame: conversion matrix from cartesian lattice [Angs-1] to current frame in [unit] below', ...
  '  rlu2frame:  conversion matrix from HKL lattice [rlu] to current frame in [unit] below', ...
  '  Q:          location of the HKL evaluation in the frame, in [unit] below', ...
  '  RM:         resolution matrix in the frame, in [unit]^3.meV', ...
  '  README:     a short description of the coordinate frame', ...
  '  unit:       the unit used for wavevector representation', ...
  '  frame:      the unit vectors defining the frame, expressed in the [frameUnit] below', ...
  '  frameUnit:  the wavevector unit used for defining the [frame]', ...
  '  frameStr:   the human readable expression of the [frame] unit vectors, in [frameUnit]', ...
  '  Bragg:      the resolution FWHM along [x,y,z,E,Vana] in [unit]^3.meV', ...
  'The available coordinate frames are:' };

for index={'rlu','cart','spec','rlu_ABC','ABC'}
  res.README{end+1} = [ sprintf('%8s', index{1}) ': ' res.(index{1}).README ];
  % complete the frameStr where empty with e.g. [frame vector][unit] for humans
  res.(index{1}).frameStr = ResLibCal_RM2RMS_Labels(res.(index{1}).frameStr, ...
                           res.(index{1}).frame, res.(index{1}).frameUnit);
end

% ------------------------------------------------------------------------------
function str = ResLibCal_RM2RMS_Labels(str, v, unit)
  % get labels for given vectors
  
  for index=1:size(v,2)
    str{index} = strtrim([ str{index} ' =[' strtrim(sprintf('%.2f ', v(:,index))) ' ' unit ']' ]);
  end

% ------------------------------------------------------------------------------
function [bragg]=rc_bragg(M)
%
% RESCAL function to calculate the widths (FWHM)
% of a Bragg peak from the resolution matrix M
%
% Called by: rc_res
% Calls  to: rc_phon
%
% Output: bragg, Qx, Qy, Qz, DEE widths, Vanadium 
% 
% ResCal5/A.T.
%
if isempty(M), bragg=[]; return; end

bragg=sqrt(8*log(2))./sqrt(diag(M));
[r,bragg(5)]=rc_phon(1,M,[0 0 0 1]); % Vanadium width: flat dispersion

bragg = bragg*2; % from hwhm to fwhm

% ------------------------------------------------------------------------------
function [rp,fwhm]=rc_phon(r0,M,C)

%
% MATLAB  routine to calculate the phonon width of a scan along a 
% vector s, and a plane defined by C.X=w. r0 is the resolution constant and
% M is the resolution matrix.
%
% A.T.
if isempty(M), rp=0; fwhm=0; return; end
T=diag(ones(4,1),0);
T(4,1:4)=C;
S=inv(T);
MP=S'*M*S;
[rp,MP]=rc_int(1,r0,MP);
[rp,MP]=rc_int(1,rp,MP);
[rp,MP]=rc_int(1,rp,MP);
fwhm=sqrt(8*log(2))/sqrt(MP(1,1));

