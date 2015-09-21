function resolution = ResLibCal_ComputeResMat(EXP)
% resolution = ResLibCal_ComputeResMat(EXP): compute the resolution function
%
% Compute the resolution matrices by calling the EXP.method
%
% Return:
%  resolution.R0:  Resolution prefactor
%  resolution.RM:  Resolution matrix with x-axis along Q
%  resolution.RMS: Resolution matrix wrt reciprocal lattice units

% Calls: ResLibCal_fig2EXP, ResMat, 
%        rc_cnmat, rc_popma, res3ax5, vTAS_AFILL, rc_bragg
%        rc_re2rc called in ResLibCal_RM2RMS

  resolution = {};
  if nargin == 0,  EXP=''; end
  if isempty(EXP), EXP=ResLibCal_fig2EXP(get(0,'CurrentFigure')); end
  if ~isstruct(EXP), return; end
  
  % check EXP structure. Perhaps it is a full ResLibCal structure
  if isfield(EXP,'EXP')
    EXP = EXP.EXP;
  end
  
  % handle case where EXP is an array or cell array
  if numel(EXP) > 1
    for index=1:numel(EXP)
      if iscell(EXP)
        resolution{index} = feval(mfilename, EXP{index});
      else
        resolution{index} = feval(mfilename, EXP(index));
      end
    end
    return
  end

  % EXP is now a single structure
  if isfield(EXP,'mono')
    if     isfield(EXP.mono,'d')   EXP.mono.tau = 2*pi/EXP.mono.d;
    elseif isfield(EXP.mono,'tau') EXP.mono.d   = 2*pi/EXP.mono.tau; end
  end
  if isfield(EXP,'ana')
    if     isfield(EXP.ana,'d')   EXP.ana.tau = 2*pi/EXP.ana.d;
    elseif isfield(EXP.ana,'tau') EXP.ana.d   = 2*pi/EXP.ana.tau; end
  end

  if isnumeric(EXP.method)
    methods = get(ResLibCal_fig('EXP_method'),'String');
    EXP.method=methods{EXP.method+1};
  end
  
  % prepare potential scan in HKLE
  QH  = EXP.QH; QK = EXP.QK; QL = EXP.QL; W =EXP.W;
  
  % handle case where all HKLE are vectors of same length
  myisvector=@(c)max(size(c)) == numel(c);
  sz = [ numel(QH) numel(QK) numel(QL) numel(W) ];
  
  if sz(1) > 1 && myisvector(QH) && myisvector(QK) && myisvector(QL) && myisvector(W) && all(sz == sz(1))
    % all are vectors of same length: line
    for index=1:sz(1)
      resolution{end+1} = ResLibCal_ComputeResMat_Single(EXP, QH(index), QK(index), QL(index), W(index));
    end
  else
    % 4D case or single scalars
  
    len = prod([ numel(QH) numel(QK) numel(QL) numel(W) ]);

    for iqh=1:numel(QH), 
    for iqk=1:numel(QK), 
    for iql=1:numel(QL), 
    for ien=1:numel(W),  
      % loop on scan steps
      res = ResLibCal_ComputeResMat_Single(EXP, QH(iqh), QK(iqk), QL(iql), W(ien));
      
      % store resolution
      if len == 1
        resolution = res;
      else
        resolution{end+1} = res;
      end
    end % for HKLE
    end
    end
    end
  end % test for vector or 4D [QH,QK,QL,W] locations
% end ResLibCal_ComputeResMat

% ------------------------------------------------------------------------------
function res= ResLibCal_ComputeResMat_Single(EXP, h,k,l,w)
% compute a single resolution at HKLE
    % initiate empty values
    R0=0; RM=[]; RMS=[]; bragg = [];

    % compute Ki, Kf, Ei, Ef
    EXP.QH = h; EXP.QK=k; EXP.QL=l; EXP.W=w;
    
    method_orig = EXP.method;
    EXP.method = lower(EXP.method);
    
    % variabe focusing for each steps
    if any([ EXP.mono.rv EXP.mono.rh EXP.ana.rv EXP.ana.rh ] <= 0)
      [rho,EXP.ki,EXP.kf] = rc_focus(EXP);  % in 'private', from ResCal5
      if EXP.mono.rv < 0, EXP.mono.rv=rho.mv; end 
      if EXP.mono.rh < 0, EXP.mono.rh=rho.mh; end 
      if EXP.ana.rv  < 0, EXP.ana.rv =rho.av; end
      if EXP.ana.rh  < 0, EXP.ana.rh =rho.ah; end
      if EXP.mono.rv==0 || ~isfinite(EXP.mono.rv), EXP.mono.rv=Inf; end
      if EXP.mono.rh==0 || ~isfinite(EXP.mono.rh), EXP.mono.rh=Inf; end
      if EXP.ana.rv==0  || ~isfinite(EXP.ana.rv),  EXP.ana.rv=Inf; end
      if EXP.ana.rh==0  || ~isfinite(EXP.ana.rh),  EXP.ana.rh=Inf; end
    end
    res.focus = [ EXP.mono.rv EXP.mono.rh EXP.ana.rv EXP.ana.rh ];
    
    % sample rotation. We ignore it in this implementation.
    % [EXP,x,y,z]       = ResLibCal_SampleRotateS(h,k,l,EXP);
    
    % choice of method
    if ~isempty(strfind(EXP.method, 'rescal5'))
      % a,b,c,alpha,beta,gamma, QH,QK,QL
      f=0.4826; % f converts from energy units into k^2, f=0.4826 for meV
      if ~isempty(strfind(EXP.method, 'cooper')), method=@rc_cnmat;
      else                                        method=@rc_popma; 
      end
      [R0,RM]   = feval(method,f,0,EXP,0);
    elseif ~isempty(strfind(EXP.method, 'res3ax'))
      if ~isempty(strfind(EXP.method, 'cooper')), method=@res3ax3;
      else                                        return; % no Popovici from J Ollivier
      end
      [R0,RM]   = feval(method,h,k,l,w, EXP);
    elseif ~isempty(strfind(EXP.method, 'afill'))
      % This method is 100% equivalent to ResCal5/Cooper-Nathans
      method    = @Rescal_AFILL; 
      [R0,RM]   = feval(method,h,k,l,w, EXP);
    else % default is 'reslib'
      method = @ResMat;
      % calls ResLib/ResMat
      [sample,rsample]=GetLattice(EXP);
      Q=modvec(h,k,l,rsample);
      [R0,RM]= ResMat(Q,w,EXP);
      % [R0, RMS, RM,x,y,z] = ResMatS(h,k,l,w, EXP);
    end
    % resolution volume and matrices -------------------------------------------
    %
    % There are 2 coordinate frames:
    % ABC is relative to the A,B vectors. These vectors should then be chosen
    %     according to the measurements that will be carried-out.
    %     The [ABC] frame is ortho-normalised (cartesian).
    % xyz is relative to Q, and vertical axis.
    %
    % RM is the resolution matrix in [xyz] frame with x // Q, z vertical
    
    res.R0    = R0;
    res.HKLE  = [ h k l w ];  % evaluation location
    res.method= method_orig;

    % resolution matrix in [abc] and transformation [HKL] -> [ABC] frame
    % S = inv([x y z]) when [x y z ] = StandardSystem(EXP);
    % S = matrix 's' in inline (below) ResLibCal_ComputeResMat_Angles
    if ~isempty(RM) && isreal(R0)
      res.spec.RM= RM;
      res                 = ResLibCal_RM2RMS(EXP, res);    % RM is other coordinate frames
      [res.angles, res.Q] = ResLibCal_ComputeResMat_Angles(h,k,l,w, ...
                                  EXP, res.ABC.rlu2frame); % spectrometer angles
      res                 = ResLibCal_RM2clouds(EXP, res); % generate MC clouds (spec, rlu)
    else
      res.README='Could not compute the resolution';
      [sample,rsample]=GetLattice(EXP);
      res.Q=modvec(h,k,l,rsample);
      disp([ datestr(now) ': ' mfilename ': ' func2str(method) ': KI,KF,Q triangle will not close (kinematic equations). ' ]);
      disp('  Change the value of KFIX,FX,QH,QK,QL or W.');
      if (res.Q < .5)
        disp('  Try an other equivalent Bragg peak further in the reciprocal lattice.');
      end
      disp([ h k l w])
    end
    
% ------------------------------------------------------------------------------
function [A,Q] = ResLibCal_ComputeResMat_Angles(h,k,l,w,EXP, s)
% compute all TAS angles (in plane)
% code extracted from TAS MAD/ILL and McStas/templateTAS.instr

    % compute angles
    fx = 2*(EXP.infin==-1)+(EXP.infin==1);
    kfix = EXP.Kfixed;
    f=0.4826; % f converts from energy units into k^2, f=0.4826 for meV
    ki=sqrt(kfix^2+(fx-1)*f*w);  % kinematical equations.
    kf=sqrt(kfix^2-(2-fx)*f*w);

    % compute the transversal Q component, and A3 (sample rotation)
    % from McStas/templateTAS.instr and TAS MAD ILL
    if nargin < 6
      a     = [ EXP.sample.a     EXP.sample.b    EXP.sample.c ]/2/pi;
      alpha = [ EXP.sample.alpha EXP.sample.beta EXP.sample.gamma ]*pi/180;
      cosa  = cos(alpha); sina = sin(alpha);
      cc    = sum(cosa.*cosa);
      cc    = 1+2*prod(cosa) - cc;
      cc    = sqrt(cc);
      b     = sina./(a*cc);
      c1    = circshift(cosa',-1); c2    = circshift(c1,-1); 
      s1    = circshift(sina',-1); s2    = circshift(s1,-1); 
      cosb  = (c1.*c2 - cosa')./(s1.*s2);
      sinb  = sqrt(1 - cosb.*cosb);

      bb    = [b(1)          0                    0 
               b(2)*cosb(3)  b(2)*sinb(3)         0
               b(3)*cosb(2) -b(3)*sinb(2)*cosa(1) 1/a(3)];
      bb = bb';
               
      aspv  = [ EXP.orient1' EXP.orient2' ];
      vv = zeros(3,3);
      vv(1:2,:)  = transpose(bb*aspv);
      for m=3:-1:2
        vt    = circshift(vv,[1 1]).*circshift(vv, [2 2]) ...
              - circshift(vv,[1 2]).*circshift(vv, [2 1]);
        vv(m,:) = vt(m,:);
      end
      c     = sqrt(sum(vv'.^2));
      vv    = vv./repmat(c,[3 1])';
      s     = vv*bb;
      % the matrix 's' is the same as 'S' from ResLibCal_RM2RMS, and its inverse 
      % is the same as [x,y,z]=StandardSystem(EXP)
    end
    
    qt    = [h k l ]*s'; % from [HKL] to [A B C] frame.
    qs    = sum(qt.*qt); Q=sqrt(qs);
    sm =EXP.mono.dir;
    ss =EXP.sample.dir;
    sa =EXP.ana.dir;
    dm =2*pi/EXP.mono.tau;
    da =2*pi/EXP.ana.tau;
    thetaa=sa*asin(pi/(da*kf));      % theta angles for analyser
    thetam=sm*asin(pi/(dm*ki));      % and monochromator.
    thetas=ss*0.5*acos((ki^2+kf^2-Q^2)/(2*ki*kf)); % scattering angle from sample.

    A3 = -atan2(qt(2),qt(1)) ...
         -ss*acos( (kf*kf-Q*Q-ki*ki)/(-2*Q*ki) );
  
    A1=thetam; A2=2*A1; A4=2*thetas; A5=thetaa; A6=2*A5;
    
    A = [A1 A2 A3 A4 A5 A6]*180/pi;

