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
    methods = get(findall(fig, 'Tag','EXP_method'),'String');
    EXP.method=methods{EXP.method+1};
  end
  
  method_orig = EXP.method;
  EXP.method = lower(EXP.method);
  
  % prepare potential scan in HKLE
  QH  = EXP.QH; QK = EXP.QK; QL = EXP.QL; W =EXP.W;
  len = max([ numel(QH) numel(QK) numel(QL) numel(W) ]);
  
  for index=1:len 
    % loop on scan steps
    if numel(QH) < index, h=QH(end); else h=QH(index); end
    if numel(QK) < index, k=QK(end); else k=QK(index); end
    if numel(QL) < index, l=QL(end); else l=QL(index); end
    if numel(W ) < index, w=W(end);  else w=W(index);  end
    % initiate empty values
    R0=1; RM=[]; RMS=[]; bragg = [];
    
    % choice of method
    if ~isempty(strfind(EXP.method, 'rescal'))
      if exist('rc_cnmat') == 2
        % a,b,c,alpha,beta,gamma, QH,QK,QL
        f=0.4826; % f converts from energy units into k^2, f=0.4826 for meV
        if ~isempty(strfind(EXP.method, 'cooper')), method=@rc_cnmat;
        else                                        method=@rc_popma; 
        end
        EXProt    = ResLibCal_SampleRotateS(h,k,l,EXP);
        EXProt.QH = h; EXProt.QK=k; EXProt.QL=l; EXProt.W=w;
        p         = ResLibCal_EXP2RescalPar(EXProt);
        [Q2c,Qmag]= rc_re2rc( p(19:21), p(22:24), p(31:33) ); 
        [R0,RM,vi,vf,Error]=feval(method,f,Qmag,p,0);
        RMS       = ResLibCal_RM2RMS(h,k,l,w,EXProt,RM);
      else
        disp([mfilename ': Rescal5/rc_cnmat is not available' ]);
      end
    elseif ~isempty(strfind(EXP.method, 'res3ax'))
      if exist('res3ax3') == 2
        if ~isempty(strfind(EXP.method, 'cooper')), method=@res3ax3;
        else                                        method=@res3ax5; 
        end
        EXProt    = ResLibCal_SampleRotateS(h,k,l,EXP);
        EXProt.QH = h; EXProt.QK=k; EXProt.QL=l; EXProt.W=w;
        [R0,RM]   = feval(method,h,k,l,w, EXProt);
        RMS       = ResLibCal_RM2RMS(h,k,l,w,EXProt,RM);
      else
        disp([mfilename ': res3ax (JO) is not available' ]);
      end
    elseif ~isempty(strfind(EXP.method, 'vtas'))
      if exist('vTAS_AFILL') == 2
        % This method is 100% equivalent to ResCal5/Cooper-Nathans
        method    = @vTAS_AFILL; 
        EXProt    = ResLibCal_SampleRotateS(h,k,l,EXP);
        EXProt.QH = h; EXProt.QK=k; EXProt.QL=l; EXProt.W=w;
        [R0,RM]   = feval(method,h,k,l,w, EXProt);
        RMS       = ResLibCal_RM2RMS(h,k,l,w,EXProt,RM);
      else
        disp([mfilename ': vTAS_AFILL is not available' ]);
      end
    else % default is 'reslib'
      if ~isempty(strfind(EXP.method, 'cooper')) EXP.method=0; 
      else                                       EXP.method=1; 
      end
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
      bragg = 2.35./sqrt(diag(RMS)); % overridden by rc_projs in rlu/Q frame when plotted
    end
    % resolution volume and matrices
    res.R0    = R0;
    res.RM    = RM;  % M in [Qx,Qy,Qz,E] frame
    res.RMS   = RMS; % M in [Q1,Q2,Qz,E] frame
    res.Bragg = bragg;
    res.HKLE  = [ h k l w ];
    res.method= method_orig;
    
    % store resolution
    if len == 1
      resolution = res;
    else
      resolution{index} = res;
    end
  end % for HKLE
% end ResLibCal_ComputeResMat

