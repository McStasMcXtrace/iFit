function res = ResLibCal_RM2RMS(EXP, res)
% RMS=ResLibCal_RM2RMS(H,K,L,W,EXP,RM) rotate Resolution matrix in [Q1 Q2] frame
% [RMS,S]=ResLibCal_RM2RMS(H,K,L,W,EXP,RM)
%   also returns the tranformation matrix in the lattice coordinates
%   Q in frame [ABC] = S*[H K L]', e.g. to pass the cloud in the ABC frame

if isempty(res.xyz.RM), return; end

% method: ResCal5 which is the same as legacy RESCAL.
% Transformation matrix 'S' checked against templateTAS.instr and TAS MAD ILL
% Results fully agree with RESCAL, McStas templateTAS and TAD MAD/ILL.

HKL = res.HKLE(1:3)';
%----- Calculate Q2c matrix
% Q2c = [a_star_vec b_star_vec c_star_vec];
[Q2c, QM, Qcart]= rc_re2rc( [ EXP.sample.a EXP.sample.b EXP.sample.c ], ...
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

%----- S transformation matrix from (h,k,l) to V1,V2,V3 = [ABC] orthonormal
S=U*Q2c;  % This is used to bring the Resolution matrix into a defined frame.

%----- Q vector in cartesian coordinates, i.e. in a ortho-normalised [a* b* c*] system.
% Qcart=Q2c*HKL;

%----- Work out angle of Q wrt to V1, V2
TT=S*HKL; % TT bring from HKL to ABC frame
cos_theta=TT(1)/norm(TT);
sin_theta=TT(2)/norm(TT);
R=[cos_theta sin_theta 0; -sin_theta cos_theta 0; 0 0 1]; % HKL in ABC frame

%----- get and transform resolution matrix in Q frame. RM is resolution matrix in Qx, Qy & Qz
T=zeros(4,4);
T(4,4)=1;
T(1:3,1:3)=R*S;

%----- Resolution ellipsoid in terms of H,K,L,EN ([Rlu] & [meV])
res.abc.RM=T'*res.xyz.RM*T;    % 'M' in ortho-normal cartesian [ABC] lattice frame = RMS...
res.abc.README = { ...
  'Frame:     [ABC] orthonormal frame expressed in [HKL rlu] as columns'
  'FrameStr:  printable representation of the [ABC] axes in HKL [rlu]'
  'RM:        resolution matrix in the [ABC] lattice frame'
  'hkl2Frame: converts [abc] = hkl2Frame*[H K L]'' coordinates'
  'Q:         location of the HKL evaluation in the [ABC] frame'
  'Bragg:     resolution FWHM along [A,B,C,E,Vana] in rlu^3.meV' };
res.abc.hkl2Frame=S;    % [Q in frame ABC ortho-normal] = hkl2Frame*[H K L]'
res.abc.Q = res.abc.hkl2Frame*HKL;       % [HKL] in the [ABC] cartesian frame



% any [HKL] vector is [Q_abc_cartesian] = (res.hkl2ABC)*[H K L]'
% res.abc.hkl2ABC      = S;    
% res.abc.hkl2ABCcart  = Q2c;

% code from e.g. ResLib/ResPlot
% the xvec,yvec,zvec computation is redundent with resolution.hkl2ABC
% [xvec,yvec,zvec,sample,rsample]=StandardSystem(EXP);

xyz  = inv(S);
xvec = xyz(:,1);  % [ABC] ortho-normal frame
yvec = xyz(:,2);
zvec = xyz(:,3);

% normalised, orthonormal ABC
o1 = EXP.orient1;
o2 = EXP.orient2;
[sample,rsample]=GetLattice(EXP);
pr = scalar(o2(1),o2(2),o2(3),yvec(1),yvec(2),yvec(3),rsample);
o2 = yvec*pr;       % o2 = yvec ortho normal

o1 = o1(:)';        % this is A
o2 = o2(:)';        % this is yvec (i.e. ortho-normal B)
o3 = cross(o1,o2);  % this is [AxB]
res.abc.FrameStr = ResLibCal_RM2RMS_Labels({o1,o2,o3});

% [ABC] ortho-normal Base vectors are columns, expressed in [HKL]
res.abc.Frame    = [ o1/norm(o1) ; o2/norm(o2) ; o3/norm(o3) ]'; 

res.abc.unit = 'rlu';

% compute 'resolution' width
res.abc.Bragg = rc_bragg(res.abc.RM); % dQ1,dQ2,dQ3,V,dE in [Q1,Q2,Q3,E] frame

% compute the reference frames -------------------------------------------------



% [xyz] ortho-normal Base vectors are columns, expressed in [HKL]
res.xyz.README = { ...
  'Frame:     [xyz] orthonormal frame expressed in [HKL rlu] as columns. x=[QH,QK,QL]'
  'FrameStr:  printable representation of the [xyz] axes in HKL [rlu]'
  'RM:        resolution matrix in the [xyz] lattice frame'
  'hkl2Frame: converts [xyz] = hkl2Frame*[H K L]'' coordinates' 
  'Q:         location of the HKL evaluation in the [xyz] frame=hkl2Frame*[H K L]'' // [1 0 0]'
  'Bragg:     resolution FWHM along [x,y,z,E,Vana] in Angs^-3.meV'};
res.xyz.hkl2Frame = T(1:3,1:3);       % [Q in frame xyz ortho-normal] = hkl2Frame*[H K L]'
res.xyz.Q = res.xyz.hkl2Frame*HKL; % by definition, should be // [1 0 0]

% code from rescal5/rc_res: give axes vectors [xyz] as HKL coordinates in lattice frame
qx=HKL;                   % this is by definition
qz=inv(Q2c)*V3;           % V3 = cross(Q2c*A,Q2c*B)
qy=cross(qz,qx);
res.xyz.FrameStr = ResLibCal_RM2RMS_Labels({qx,qy,qz});

qx=qx/norm(qx);
qy=qy/norm(qy);
qz=qz/norm(qz);
res.xyz.Frame    = [ qx qy qz ]; % xyz vectors in HKL frame, as columns

res.xyz.unit = 'Angs-1';
res.xyz.Bragg  = rc_bragg(res.xyz.RM); % dQx,dQy,dQz,V,dE in [Qx,Qy,Qz,E] frame


% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% method2: ResLib
% the code extracted from ResLib/ResMatS does not seem to properly compute the 
% RMS matrix, as it assumes scattering to take place in plane, and restricts
% the transformation matrix to [1:2] sub matrix.
% Beware: the code assumes also rows-columns 3 and 4 have been swapped in RM.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
% ------------------------------------------------------------------------------
function str = ResLibCal_RM2RMS_Labels(v)
  % get labels for given vectors
  
  for index=1:numel(v)
    o2 = v{index}; o2=o2(:)'; ro=round(o2); so=round(sum(o2.*o2));
    % test if the normalised vector has a simple rational expression
    if all(abs(o2 - ro) < 1e-2)
      o2 = ro;
      if so > 1,  o2=[ '[' strtrim(sprintf('%i ', o2)) ']/\surd' num2str(so) ];
      else        o2=[ '[' strtrim(sprintf('%i ', o2)) ']' ];
      end
    else o2 = [ '[' strtrim(sprintf('%.2f ', o2)) ']' ];
    end
    str{index} = o2;
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

