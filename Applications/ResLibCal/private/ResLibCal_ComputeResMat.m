function resolution = ResLibCal_ComputeResMat(EXP)
% resolution = ResLibCal_ComputeResMat(EXP): compute the resolution function
%
% Compute the resolution matrices by calling the EXP.method
%
% Return:
%  resolution.R0:  Resolution prefactor
%  resolution.RM:  Resolution matrix with x-axis along Q
%  resolution.RMS: Resolution matrix wrt reciprocal lattice units

% Calls: ResLibCal_fig2EXP, ResMatS, ResLibCal_SampleRotateS, rc_re2rc
%        rc_cnmat, rc_popma, res3ax5, vTAS_AFILL, rc_bragg

  resolution = [];
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
  
  method_orig = EXP.method;
  EXP.method = lower(EXP.method);
  
  % prepare potential scan in HKLE
  QH  = EXP.QH; QK = EXP.QK; QL = EXP.QL; W =EXP.W;
  len = prod([ numel(QH) numel(QK) numel(QL) numel(W) ]);
  
  EXPorg = EXP;

  for iqh=1:numel(QH), h = QH(iqh);
  for iqk=1:numel(QK), k = QK(iqk);
  for iql=1:numel(QL), l = QL(iql);
  for ien=1:numel(W),  w = W(ien);
    % loop on scan steps
    
    % initiate empty values
    R0=1; RM=[]; RMS=[]; bragg = [];

    % compute Ki, Kf, Ei, Ef
    EXP=EXPorg;
    EXP.QH = h; EXP.QK=k; EXP.QL=l; EXP.W=w;
    
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
    res.RMV=EXP.mono.rv;
    res.RMH=EXP.mono.rh;
    res.RAV=EXP.ana.rv;
    res.RAH=EXP.ana.rh;
    
    % choice of method
    if ~isempty(strfind(EXP.method, 'rescal5'))
      if exist('rc_cnmat') == 2
        % a,b,c,alpha,beta,gamma, QH,QK,QL
        f=0.4826; % f converts from energy units into k^2, f=0.4826 for meV
        if ~isempty(strfind(EXP.method, 'cooper')), method=@rc_cnmat;
        else                                        method=@rc_popma; 
        end
        EXP    = ResLibCal_SampleRotateS(h,k,l,EXP);
        % p         = ResLibCal_EXP2RescalPar(EXP);
        % [Q2c,Qmag]= rc_re2rc( p(19:21), p(22:24), p(31:33) ); 
        [R0,RM]   = feval(method,f,0,EXP,0);
        RMS       = ResLibCal_RM2RMS(h,k,l,w,EXP,RM);
      else
        disp([mfilename ': Rescal5/rc_cnmat is not available' ]);
      end
    elseif ~isempty(strfind(EXP.method, 'res3ax'))
      if exist('res3ax3') == 2
        if ~isempty(strfind(EXP.method, 'cooper')), method=@res3ax3;
        else                                        return; % no Popovici from J Ollivier
        end
        EXP    = ResLibCal_SampleRotateS(h,k,l,EXP);
        [R0,RM]   = feval(method,h,k,l,w, EXP);
        RMS       = ResLibCal_RM2RMS(h,k,l,w,EXP,RM);
      else
        disp([mfilename ': res3ax (JO) is not available' ]);
      end
    elseif ~isempty(strfind(EXP.method, 'afill'))
      if exist('Rescal_AFILL') == 2
        % This method is 100% equivalent to ResCal5/Cooper-Nathans
        method    = @Rescal_AFILL; 
        EXP    = ResLibCal_SampleRotateS(h,k,l,EXP);
        [R0,RM]   = feval(method,h,k,l,w, EXP);
        RMS       = ResLibCal_RM2RMS(h,k,l,w,EXP,RM);
      else
        disp([mfilename ': rescal/AFILL is not available' ]);
      end
    else % default is 'reslib'
      % calls ResLib/ResMatS
      % depends: ResMatS, CleanArgs, StandardSystem
      %          modvec, scalar, ResMat, GetLattice, star, GetTau
      if exist('ResMatS') == 2
        [R0, RMS, RM] = ResMatS(h,k,l,w, EXP);
      else
        disp([mfilename ': ResLib 3.4/ResMatS is not available' ]);
      end
    end
    % assemble 'resolution' step
    if ~all(isreal(RM))  RM=[]; end
    if ~all(isreal(RMS)) RMS=[]; end
    if ~isempty(RMS)
      % compute some widths in [Q1,Q2,Qz,E]
      res.BraggS = rc_bragg(RMS); % dQ1,dQ2,dQz,V,dE in [Q1,Q2,Qz,E] frame
    end
    if ~isempty(RM)
      % compute some widths in [Qx,Qy,Qz,E]
      res.Bragg = rc_bragg(RM); % dQx,dQy,dQz,V,dE in [Qx,Qy,Qz,E] frame
    end
    % resolution volume and matrices
    res.R0    = R0;
    res.RM    = RM;  % M in [Qx,Qy,Qz,E] frame
    res.RMS   = RMS; % M in [Q1,Q2,Qz,E] frame
    res.HKLE  = [ h k l w ];
    res.method= method_orig;
    [res.angles, res.Q]     = ResLibCal_ComputeResMat_Angles(h,k,l,w,EXP);
    
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
% end ResLibCal_ComputeResMat

% ------------------------------------------------------------------------------
function [A,Q] = ResLibCal_ComputeResMat_Angles(h,k,l,w,EXP)
% compute all TAS angles (in plane)

    % compute angles
    fx = 2*(EXP.infin==-1)+(EXP.infin==1);
    kfix = EXP.Kfixed;
    f=0.4826; % f converts from energy units into k^2, f=0.4826 for meV
    ki=sqrt(kfix^2+(fx-1)*f*w);  % kinematical equations.
    kf=sqrt(kfix^2-(2-fx)*f*w);

    % compute the transversal Q component, and A3 (sample rotation)
    % from McStas templateTAS.instr and TAS MAD ILL
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
    qt    = [h k l ]*s';
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

