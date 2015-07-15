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
    res.focus = [ EXP.mono.rv EXP.mono.rh EXP.ana.rv EXP.ana.rh ];
    
    % sample rotation. We ignore it in this implementation.
    % [EXP,x,y,z]       = ResLibCal_SampleRotateS(h,k,l,EXP);
    
    % choice of method
    if ~isempty(strfind(EXP.method, 'rescal5'))
      if exist('rc_cnmat') == 2
        % a,b,c,alpha,beta,gamma, QH,QK,QL
        f=0.4826; % f converts from energy units into k^2, f=0.4826 for meV
        if ~isempty(strfind(EXP.method, 'cooper')), method=@rc_cnmat;
        else                                        method=@rc_popma; 
        end
        [R0,RM]   = feval(method,f,0,EXP,0);
      else
        disp([mfilename ': Rescal5/rc_cnmat is not available' ]);
      end
    elseif ~isempty(strfind(EXP.method, 'res3ax'))
      if exist('res3ax3') == 2
        if ~isempty(strfind(EXP.method, 'cooper')), method=@res3ax3;
        else                                        return; % no Popovici from J Ollivier
        end
        [R0,RM]   = feval(method,h,k,l,w, EXP);
      else
        disp([mfilename ': res3ax (JO) is not available' ]);
      end
    elseif ~isempty(strfind(EXP.method, 'afill'))
      if exist('Rescal_AFILL') == 2
        % This method is 100% equivalent to ResCal5/Cooper-Nathans
        method    = @Rescal_AFILL; 
        [R0,RM]   = feval(method,h,k,l,w, EXP);
      else
        disp([mfilename ': rescal/AFILL is not available' ]);
      end
    else % default is 'reslib'
      % calls ResLib/ResMat
      if exist('ResMat') == 2
        [sample,rsample]=GetLattice(EXP);
        Q=modvec(h,k,l,rsample);
        [R0,RM]= ResMat(Q,w,EXP);
        % [R0, RMS, RM,x,y,z] = ResMatS(h,k,l,w, EXP);
        
      else
        disp([mfilename ': ResLib 3.4/ResMat is not available' ]);
      end
    end
    % resolution matrix in rlu and transformation [HKL] -> [ABC] frame
    % S = inv([x y z]) when [x y z ] = StandardSystem(EXP);
    % S = matrix 's' in inline (below) ResLibCal_ComputeResMat_Angles
    [RMS,S, U]   = ResLibCal_RM2RMS(h,k,l,w,EXP,RM);
    
    % assemble 'resolution' step
    if ~all(isreal(RM))  RM=[]; end
    if ~all(isreal(RMS)) RMS=[]; end
    if ~isempty(RMS)
      % compute some widths in [Q1,Q2,Q3,E]
      res.BraggS = rc_bragg(RMS); % dQ1,dQ2,dQz,V,dE in [Q1,Q2,Q3,E] frame
    end
    if ~isempty(RM)
      % compute some widths in [Qx,Qy,Qz,E]
      res.Bragg = rc_bragg(RM); % dQx,dQy,dQz,V,dE in [Qx,Qy,Qz,E] frame
    end
    % resolution volume and matrices
    res.R0    = R0;
    res.RM    = RM;  % M in [Qx,Qy,Qz,E] frame Qx // Q
    res.RMS   = RMS; % M in [QA,QB,QC,E] frame
    res.HKLE  = [ h k l w ];
    % res.vectors=[ x y z ];
    % vectors = inv(hkl2ABC)
    res.hkl2ABC= S; % [hkl in ABC frame] = (res.hkl2ABC)*[ h k l ]'
    res.U      = U; % ortho-normal ABC frame
    res.method = method_orig;
    [res.angles, res.Q]     = ResLibCal_ComputeResMat_Angles(h,k,l,w,EXP, S);
    
    [EXP, res] = ResLibCal_ComputeResMat_Frame(EXP, res);
    
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

% ------------------------------------------------------------------------------
function [EXP, resolution] = ResLibCal_ComputeResMat_Frame(EXP, resolution)
  
  % the xvec,yvec,zvec computation is redundent with resolution.hkl2ABC
  % [xvec,yvec,zvec,sample,rsample]=StandardSystem(EXP);

  [sample,rsample]=GetLattice(EXP);
  xyz = inv(resolution.hkl2ABC);
  xvec = xyz(:,1);
  yvec = xyz(:,2);
  zvec = xyz(:,3);
  
  % normalised, orthonormal ABC
  o1=EXP.orient1;
  o2=EXP.orient2;
  pr=scalar(o2(1),o2(2),o2(3),yvec(1),yvec(2),yvec(3),rsample);
  o2 = yvec*pr; % o2 // yvec ortho normal

  o1 = o1(:)';
  o2 = o2(:)';
  o3 = cross(o1,o2);

  resolution.rluFrame = [ o1 ; o2 ; o3 ]'; % Base vectors are columns
  
  % now get labels
  % convert o1 and o2 to normalised strings
  if all(abs(o1 - round(o1)) < 1e-5)
    o1 = round(o1);
    if sum(o1.*o1) > 1,  o1=[ '[' sprintf('%i ', o1) ']/\surd' num2str(sum(o1.*o1)) ];
    else                 o1=[ '[' sprintf('%i ', o1) ']' ];
    end
  else o1 = [ '[' sprintf('%.2f ', o1) ']' ]; end

  if all(abs(o2 - round(o2)) < 1e-2)
    o2 = round(o2);
    if sum(o2.*o2) > 1,  o2=[ '[' sprintf('%i ', o2) ']/\surd' num2str(sum(o2.*o2)) ];
    else                 o2=[ '[' sprintf('%i ', o2) ']' ];
    end
  else o2 = [ '[' sprintf('%.2f ', o2) ']' ]; end
  
  if all(abs(o3 - round(o3)) < 1e-5)
    o3 = round(o3);
    if sum(o3.*o3) > 1,  o3=[ '[' sprintf('%i ', o3) ']/\surd' num2str(sum(o3.*o3)) ];
    else                 o3=[ '[' sprintf('%i ', o3) ']' ];
    end
  else o3 = [ '[' sprintf('%.2f ', o3) ']' ]; end
  
  resolution.rluFrameStr = { o1, o2, o3 };

