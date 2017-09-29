function [resolution,R0,RM] = ResLibCal_RM2clouds_mcstas(EXP, resolution, angles)
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

R0=0; RM=[];

if ischar(EXP) && strcmp(EXP,'compile')
  resolution = ResLibCal_compile_mcstas_tas;
  return
end

% accuracy obtrained from McStas estimates:
% NMC=1000  -> 10%
% NMC=10000 -> 2.5%
NMC=get(ResLibCal_fig('View_NMC'), 'UserData');
if     isfield(EXP, 'NMC'),         NMC=EXP.NMC;
elseif isfield(resolution, 'NMC'),  NMC=resolution.NMC; end
if isempty(NMC), NMC  = 200; end

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
% WARNING: angles are given with signs, so we use abs(angles) and
%'SM',p.SM, 'SS',p.SS, 'SA',p.SA, ... -> 1
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
  'radius',2*p.WS/100, 'height',p.HS/100, ...
  'WB',p.WB/100,'HB',p.HB/100, ...
  'WD',p.WB/100,'HD',p.HB/100, ...
  'A1',abs(angles(1)),'A2',abs(angles(2)), ...
  'A3',    angles(3), 'A4',abs(angles(4)), ...
  'A5',abs(angles(5)),'A6',abs(angles(6)));

% now launch the Mcstas simulation with given parameters
if ismac,      precmd = 'DYLD_LIBRARY_PATH= ;';
elseif isunix, precmd = 'LD_LIBRARY_PATH= ; '; 
else precmd=''; end

% get temporary directory for results
d = tempname;
cmd = sprintf(' --dir=%s -n %i ',d, NMC*400);

for f=fieldnames(pars)'
  cmd = [ cmd ' ' f{1} '=' num2str(pars.(f{1})) ];
end
disp([ compiled cmd ]);
[status, result] = system([ precmd compiled cmd ]);

% read results
if ~isdir(d)
  fprintf(1, '%s: Failed running McStas TAS model %s\n', mfilename, compiled);
  return
end

if ~isempty(dir(fullfile(d,'resolution.dat')))
  % read the cloud in Angs-1, independent of the sample lattice, except A3
  cloud = iLoad(fullfile(d,'resolution.dat'),'mccode');
  p  = prod(cloud.data(:,10:11),2); % pi*pf = intensity (e.g. gaussian profile)
  R0 = sum(p)*1e6;  % total intensity
  
  % we select only points which are within the RMS. The central part (p > 50%) is 
  % always extracted, the lower part follows a random distribution.
  threshold = max(p)*rand(size(p))/2;
  index     = find(p >= threshold);
  
  Ki= cloud.data(index,1:3); 
  Kf= cloud.data(index,4:6);
  Q = Ki-Kf;  % Q transfer in the A3=<Ki,A> rotated frame, e.g. ABC lattice (Ang-1) frame
  % QM= sqrt(sum(Q.^2,2));  % not needed
  VS2E= 5.22703725e-6;  % from McStas
  K2V = 629.622368;
  Vi = K2V*sqrt(sum(Ki.^2,2));
  Vf = K2V*sqrt(sum(Kf.^2,2));
  Ei = VS2E*Vi.^2;
  Ef = VS2E*Vf.^2;
  E  = Ei-Ef; % energy transfer
  
  % McStas orientation
  % Z along A3, X vertical, Y on the 'right' from Z (looking at analyser).
  % ResLibCal   McStas
  % X           Z
  % Y           X
  % Z           Y
  HKLE = Q(:,[3 1 2])';          % in [ABC]
  resolution.ABC.cloud = { HKLE(1,:)' HKLE(2,:)' HKLE(3,:)' E };
  HKLE_ABC = HKLE;

  % compute the cloud for other frames: 
  
  % [ABC] U -> [rlu] R
  U2R = inv(resolution.ABC.rlu2frame);
  HKLE = U2R*HKLE;  % in [rlu] from [ABC] by ABC.frame2rlu
  resolution.rlu.cloud = { HKLE(1,:)' HKLE(2,:)' HKLE(3,:)' E };

  % [ABC] U -> [spec] Q
  % U2Q = R2Q * U2R = res.spec.rlu2frame * inv(res.ABC.rlu2frame)
  U2Q = resolution.spec.rlu2frame * inv(resolution.ABC.rlu2frame);
  HKLE = U2Q*HKLE_ABC;      % in [spec] from [ABC]
  resolution.spec.cloud= { HKLE(1,:)' HKLE(2,:)' HKLE(3,:)' E };
  
  % [ABC] U -> [a*b*c* cartesian] B
  % U2B = U;  % ortho-normal res.ABC.cart2frame
  U = resolution.ABC.cart2frame';
  U2B  = U;
  HKLE = U2B*HKLE_ABC;      % in [cart] from [ABC]
  resolution.cart.cloud= { HKLE(1,:)' HKLE(2,:)' HKLE(3,:)' E };
  
  disp([ mfilename ': using ' num2str(numel(E)) ' points in cloud out of ' num2str(numel(p)) ]);
  p = p(index);

  % compute the resolution matrix on the central part (FWHM), starting from [ABC]
  % index  = find(p > max(p)*.65);
  HKLE   = [ resolution.ABC.cloud{:} ];
  % there is an issue with the apparent ellipsoids when we switch randomly from 
  % the 'try' (few points) or 'catch' (all points). So now we use the full set.
  RM_U = MinVolEllipse(HKLE'); % in [ABC]
  
  % now we convert to [rlu]
  % RM_U = T'*RM_Q*T
  % T=diag([0 0 0 1]); T(1:3,1:3)=Q'*U;
  % res.ABC.cart2frame    = U'
  % res.spec.cart2frame   = Q'
  Q = resolution.spec.cart2frame';
  T=diag([0 0 0 1]); T(1:3,1:3)=Q'*U;
  RM = T*RM_U*T'; % in the cartesian spectrometer space
else
  fprintf(1,'%s: Failed running McStas TAS model %s: no neutron counts on detector (R0=0)\n', mfilename, compiled);
end
rmdir(d, 's');

% ------------------------------------------------------------------------------
% -------------------------------------------------------------------------
function compiled = ResLibCal_compile_mcstas_tas(compile)
  % compile templateTAS as binary when does not exist yet
  
  compiled = ''; 
  if nargin == 0, compile = ''; end
  if ismac,      precmd = 'DYLD_LIBRARY_PATH= ;';
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
  
  % when we get there, compile templateTAS_arch, not existing yet
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
      if strcmp(cc,'mpicc'), cmd{end+1} = '-DUSE_MPI'; end
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
  

