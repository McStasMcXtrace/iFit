function resolution = ResLibCal_RM2clouds_mcstas(EXP, resolution)
% ResLibCal_RM2clouds_mcstas: creates an axis system of Monte-Carlo points
%   which represent the resolution function in [abc] and [xyz] frames, using 
%   McStas templateTAS
%
%   EXP: structure containing application configuration. When not given, extracts
%        it from the main window.
%   resolution: the resolution computation structure, for both [abc] and [xyz] frames
%
% Returns:
%   resolution.rlu.cloud:  cloud of points 4D axes as { H,K,L,W } in [ABC] frame
%   resolution.spec.cloud: cloud of points 4D axes as { H,K,L,W } in [xyz] frame

% insprired from ResCal5/mc_conv
% this routine uses: in [abc|xyz] frame: RM (that is along [ABC|xyz] axes)
%                    the HKLE location, converted into [ABC|xyz] frame with hkl2Frame)

persistent labels_c;
persistent compiled;

% accuracy obtrained from McStas estimates:
% NMC=1000  -> 10%
% NMC=10000 -> 2.5%
NMC=get(ResLibCal_fig('View_NMC'), 'UserData');
if     isfield(EXP, 'NMC'),         NMC=EXP.NMC;
elseif isfield(resolution, 'NMC'),  NMC=resolution.NMC; end
if isempty(NMC), NMC  = 2000; end

% get Rescal parameters
[p, labels] = ResLibCal_EXP2RescalPar(EXP);
if isempty(compiled) % only the first time it starts or when explicitly requested
  compiled = ResLibCal_compile_mcstas_tas;  % check if templateTAS is there, or compile it
end

if ~isempty(p)
  % update equivalent RESCAL parameters from EXP
  p = mat2cell(p(:),ones(1,length(p)));
  if isempty(labels_c), labels_c=strtok(labels); end
  p = cell2struct(p(:),labels_c(:),1);
end

% assemble the McStas command line arguments as a struct
pars=struct( ...
  'KI',EXP.ki, ...
  'EN',EXP.W, ...
  'L1',p.L1/100, 'L2',p.L2/100, 'L3',p.L3/100, 'L4',p.L4/100, ...
  'SM',p.SM, 'SS',p.SS, 'SA',p.SA, ...
  'DM',p.DM, 'DA',p.DA, ...
  'RMV',1/p.ROMV, 'RMH',1/p.ROMH, 'RAV',1/p.ROAV, 'RAH',1/p.ROAH, ...
  'ETAM',p.ETAM, 'ETAA',p.ETAA, ...
  'ALF1',p.ALF1, 'ALF2',p.ALF2, 'ALF3',p.ALF3, 'ALF4',p.ALF4, ...
  'BET1',p.BET1, 'BET2',p.BET2, 'BET3',p.BET3, 'BET4',p.BET4, ...
  'NHM',9, 'NVM',9, 'WM',p.WM/100, 'HM',p.HM/100, ...
  'NHA',9, 'NVA',9, 'WA',p.WA/100, 'HA',p.HA/100, ...
  'radius',2*p.WS, 'height',p.HS, ...
  'A1',resolution.angles(1),'A2',resolution.angles(2), ...
  'A3',resolution.angles(3),'A4',resolution.angles(4), ...
  'A5',resolution.angles(5),'A6',resolution.angles(6));

% now launch the Mcstas simulation with given parameters
if ismac,  precmd = 'DYLD_LIBRARY_PATH= ;';
elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
else precmd=''; end

cmd = '';
for f=fieldnames(pars)'
  cmd = [ cmd ' ' f{1} '=' num2str(pars.(f{1})) ];
end
disp([ compiled cmd ]);
[status, result] = system([ precmd compiled cmd ]);


%----- 

% method: rescal5/rc_conv
% this code is very compact and efficient, after re-factoring and testing 
% against rescal5.

RMC = randn(4,NMC); % Monte Carlo points

% [rlu] [R] frame: Resolution ellipsoid in terms of H,K,L,EN ([Rlu] & [meV])
for frames={'rlu','spec','ABC'}  % others: 'cart','rlu_ABC','ABC'
  frame = resolution.(frames{1});
  M=frame.RM;
  [V,E]=eig(M);
  sigma=1./sqrt(diag(E)); % length along principal axes of gaussian

  % compute MC points on axes (with NMC points)
  xp   = bsxfun(@times,sigma,RMC);
  % get cloud in HKLE [rlu^3.meV], centred at 0.
  XMC  = inv(V)'*xp; % this is delta(HKLE) as a Gaussian distribution
  % compute the HKLE position in the lattice frame [a*,b*,c*,w]
  HKLE(1:3) = frame.rlu2frame*resolution.HKLE(1:3)';
  HKLE(4)   = resolution.HKLE(4);
  HKLE = bsxfun(@plus,HKLE',XMC); % add HKLE location to Gaussian

  resolution.(frames{1}).cloud = { HKLE(1,:)' HKLE(2,:)' HKLE(3,:)' HKLE(4,:)' }; % get 1D arrays per axis
  clear HKLE
end

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% method: ResLib/ConvRes
% the code extracted from ConvRes does not seem to produce sensible Monte-Carlo
% points. The distribution is clearly not Gaussian, and extends very far from 
% the HKLE position.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



% and now we can evaluate the function onto the axes.... and sum all values
% sum(feval(model, parameters, ax{:}))*resolution.R0/NMC

% the opposite operation (cloud -> RM) is computed from:
% <http://stackoverflow.com/questions/3417028/ellipse-around-the-data-in-matlab>
% B as columns
% Center = mean(B,1);
% X0 = bsxfun(@minus,B,Center);
% RM=X0'*X0 ./ (size(X0,1)-1);

% ------------------------------------------------------------------------------
% -------------------------------------------------------------------------
function compiled = ResLibCal_compile_mcstas_tas(compile)
  % compile looktxt as binary when does not exist yet
  
  compiled = ''; 
  if nargin == 0, compile = ''; end
  if ismac,  precmd = 'DYLD_LIBRARY_PATH= ;';
  elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
  else precmd=''; end
  
  if ispc, ext='.exe'; else ext=''; end
  this_path = fileparts(which(mfilename));
  
  % try in order: global(system), local, local_arch
  for try_target={ ...
          fullfile(this_path, [ 'templateTAS_' computer('arch') ext ]), ...
          fullfile(this_path, [ 'templateTAS' ext ]), ...
          [ 'templateTAS' ext ], 'templateTAS'}
      
    [status, result] = system([ try_target{1} ' --help' ]); % run from Matlab

    if status == 0 && nargin == 0
        % the executable is already there. No need to make it .
        compiled = try_target{1};
        return
    end
  end
  
  % when we get there, compile looktxt_arch, not existing yet
  target = fullfile(this_path, [ 'templateTAS_' computer('arch') ext ]);

  % search for a C compiler
  cc = '';
  for try_cc={getenv('CC'),'cc','gcc','ifc','pgcc','clang','tcc'}
    if ~isempty(try_cc{1})
      [status, result] = system([ precmd try_cc{1} ]);
      if status == 4 || ~isempty(strfind(result,'no input file'))
        cc = try_cc{1};
        break;
      end
    end
  end
  if isempty(cc)
    if ~ispc
      disp([ mfilename ': ERROR: C compiler is not available from PATH:' ])
      disp(getenv('PATH'))
      disp([ mfilename ': You may have to extend the PATH with e.g.' ])
      disp('setenv(''PATH'', [getenv(''PATH'') '':/usr/local/bin'' '':/usr/bin'' '':/usr/share/bin'' ]);');
    end
    error('%s: Can''t find a valid C compiler. Install any of: gcc, ifc, pgcc, clang, tcc\n', ...
    mfilename);
  else
    try
      fprintf(1, '%s: compiling templateTAS binary (using %s)...\n', mfilename, cc);
      cmd={cc, '-O2','-o',target, ...
        fullfile(this_path,'templateTAS.c'),'-lm'};
      cmd = sprintf('%s ',cmd{:});
      disp(cmd)
      [status, result] = system([ precmd cmd ]);
      if status == 0
        compiled = target;
      end
    end
  end

  if isempty(compiled) && ~isempty(compile)
    error('%s: Can''t compile templateTAS.c binary\n       in %s\n', ...
        mfilename, fullfile(this_path));
  end
  

