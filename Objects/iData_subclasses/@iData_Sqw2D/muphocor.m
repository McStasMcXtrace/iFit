function g = muphocor(s, varargin)
  % run MUPHOCOR on the given S(q,w) data set, incident wavelength and Temperature
  %
  %   g = muphocor(s, lambda, T, amasi, sigi, conci)
  %
  % When not given, the physical parameters are searched in the object.
  % The material physical properties include the masses, the incoherent 
  %   scattering cross sections, and the concentrations which can be given as 
  %   vectors per element/non equivalent atom in the material.
  % The concentrations are the proportion of atoms in the material.
  %
  % For H2O, one should specify manually:
  %   'amasi',  [  1.0079     15.9994 ]
  %   'sigi',   [ 82.0168      4.2325 ]
  %   'conci',  [  2           1      ]
  %
  %  Input arguments can be given in order, or with name-value pairs, or as a 
  %    structure with named fields.
  %
  %
  % MUPHOCOR Parameters:
  % --------------------
  %
  % Can be given as name/value pairs.
  %  
  %      E0:        INC. ENERGY(MEV)             
  %      FP:        FLIGHTPATH(CM)
  %      XNEL:      CHANNEL OF ELASTIC LINE    
  %      CW:        CHANNEL WIDTH(MYS)
  %      FIMI:      MIN.SCATT.ANGLE            
  %      FIMA:      MAX.SCATT.ANGLE
  %      UNT:       CONST.BACKGROUND            
  %      ABK:       Absorption coefficient
  %      AEMP:      Coeff. for calculating the detector efficiency bac=unt/(1.0-EXP(-aemp/SQRT(e0-hot(n))))
  %      NSPEC:       NUMBER OF DIFFERENT ATOMIC SPECIES
  %      HOX:       UPPER LIMIT OF DENSITY OF STATES
  %      TEMPO:     TEMPERATURE IN KELVIN
  %      DW:        MEAN DEBY WALLER COEFFICIENT
  %      AMASI:     ATOMIC MASS (g/mol)
  %      SIGI :     SIGMA (total scattering cross section, barns)
  %      CONCI:     CONCENTRATION
  %      ALFI:      SCATTERING POWER
  %      NPHO:      NUMBER OF MULTI PHONON TERMS
  %      ITM:       TOTAL NUMBER OF ITERATIONS
  %      IVIT=0(1): ITERATION BY DIFFERENCE(QUOTIENT) METHOD
  %      IDW=0:     DW COEFF. KEPT CONST.=INPUT VALUE,DW=1:DW COEFF.ITERATED
  %      IEMP=1(0): CORRECTIONS FOR COUNTER EFFICIENCY WILL BE(NOT BE) DONE
  %      IRES=1(0): DATA WILL BE(NOT BE) CORRECTED FOR SPECTR. RESOLUTION
  %      IPR=1(0):  Convolutions integrals o the multiphonon term are printed
  %      ILOSS=1:   Analysis OF the imcomplete energy loss spectrum will be done
  %      NU:        First channel OF spektrum NO: Last channel OF spektrum
  %      NUU:       First channel OF TOF distribution used FOR calculation
  %      NOO:       Last channel  OF TOF distribution used FOR calculation
  %      IGLU:      Number OF smoothing processes in CASE OF a time-dependant background
  %      UNT:       CONSTANT BACKGROUND
  %      FUN:       MULTIPLICATION FACTOR FOR TIME DEPENDENT BACKGROUND
  %      Z00:       TIME OF FLIGHT DISTRIBUTION
  %      UN0:       TIME DEPENDENT BACKGROUND (if fun = 1)
  %
  % input:
  %   s:      iData_Sqw2D object S(q,w)
  %   lambda:  optional incident neutron wavelength [Angs]
  %   temp0:   optional Temperature [K]
  %   amasi:   optional material mass [g/mol]
  %   sigi:    optional total scattering cross section [barns]
  %   conci:   optional material concentration [1]
  %
  % References:
  %   H. Schober, Journal of Neutron Research 17 (2014) 109–357
  %     DOI 10.3233/JNR-140016 (see esp. pages 328-331)
  %   V.S. Oskotskii, Sov. Phys. Solid State 9 (1967), 420.
  %   A. Sjolander, Arkiv for Fysik 14 (1958), 315.
  %   W. Reichardt, MUPHOCOR Karlsruhe Report 13.03.01p06L (1984)
  
  
  % check if muphocor is compiled, else compile it
  persistent compiled
  
  if isempty(compiled)
    compiled = muphocor_compile_binary; % check and possibly compile MUPHOCOR. Return path the EXE
  end
  
  % store parameters for iData_Sqw2D => MUPHOCOR. Lower case parameters from varargin
  p = varargin2struct({'lambda' 'temp0' 'amasi' 'sigi' 'conci'}, varargin, true);
  
  if isempty(p.temp0) && isfield(p, 't')
    p.temp0 = p.t;  % alias for temperature
  end
  if isempty(p.lambda) && isfield(p, 'e0')
    p.lambda = sqrt(81.805./p.e0);
  end
  
  if isempty(p.lambda)
    [s,p.lambda,p.distance,p.cw,~,p.wavevector] = Sqw_search_lambda(s);
  else
    [s,~,p.distance,p.cw,~,p.wavevector]      = Sqw_search_lambda(s);
  end
  if ~isfield(p, 'e0')
    p.e0    = 81.805./p.lambda.^2;
  end
  if p.cw < 1e-4
    p.cw = p.cw*1e6;  % [s] to [us];
  end
  
  if isempty(p.temp0)
    p.temp0    = Sqw_getT(s);
  end
  if isempty(p.temp0)
    p.temp0 = 293;
    disp([ mfilename ': WARNING: Using Temperature=' num2str(p.temp0) ' [K] for data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
  end
  
  % check for other required parameters: mass/we0ght, b_inc, concentration
  if isempty(p.amasi)
    p.amasi = Sqw_getT(s, {'weight'	'mass' 'AWR' 'amasi'}, 'raw');
  end
  if isempty(p.amasi)
    disp( [ mfilename ': Set object ' inputname(1) '.weight to the total molar weight' ])
    disp( [ '  or a vector with molar weight per atom in the material' ]);
    disp( [ '  or an "amasi" (4th) input argument to ' mfilename ]);
    error([ mfilename ': ERROR: Unspecified weight [g/mol] in data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
  end
  if isempty(p.sigi)
    inc = Sqw_getT(s, {'sigma_inc'}, 'raw');
    coh = Sqw_getT(s, {'sigma_coh', 'sigma', 'sigi'}, 'raw');
    p.sigi = inc+coh; % total scattering cross section
  end
  if isempty(p.sigi)
    disp( [ mfilename ': Set object ' inputname(1) '.sigma_inc to the total incoherent neutron scattering cross section' ])
    disp( [ '  or a vector with sigma_inc per atom in the material' ]);
    disp( [ '  or a "sigi" (5th) input argument to ' mfilename ]);
    error([ mfilename ': ERROR: Unspecified sigma_inc [barn] in data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
  end
  if isempty(p.conci)
    p.conci = Sqw_getT(s, {'concentration'}, 'raw');
  end
  if isempty(p.conci) || numel(p.conci) ~= numel(p.amasi)
    p.conci = ones(size(p.amasi));
  end
  
  % MUPHOCOR only works with 1 or 2 species (nspec), and then requires the partials.
  % we use it here in integrated mode.
  p.conci = p.conci/sum(p.conci);
  p.sigi  = sum(p.conci .* p.sigi);
  alpha   = p.conci .* p.sigi ./ p.amasi;
  alpha   = alpha / sum(alpha);
  p.amasi = sum(alpha .* p.amasi);  % Reichardt eq (2.13-2.15)

  % test if data set is only w > 0 or w < 0. Then print a warning and symmetrize, and Bosify
  w    = getaxis(s, 1);
  q    = getaxis(s, 2);
  if all(w(:) >=0) || all(w(:) <= 0)
    disp([ mfilename ': WARNING: the data set is given only on an energy side. Symmetrizing and applying Temperature Bose factor.' ])
    s = symmetrize(s);
    s = Bosify(s, p.temp0);
    w = getaxis(s, 1);
  end
  
  % data must be [\int sin(theta) S(theta,t) dtheta] that is start from time distribution
  % but the computation in q is much faster as we enter a (q,w) data set.
  SE2V = 437.393377;        % Convert sqrt(E)[meV] to v[m/s]
  V2K  = 1.58825361e-3;     % Convert v[m/s] to k[1/AA]
  K2V  = 1/V2K;
  VS2E = 5.22703725e-6;     % Convert (v[m/s])**2 to E[meV]
  Ki   = 2*pi/p.lambda;
  Vi   = K2V*Ki;
  Ei   = VS2E*Vi.^2;
  
  Ef   = Ei - w;
  Vf   = sqrt(Ef/VS2E); Vf(Ef<=0) = nan;
  Kf   = Vf/K2V;
  St   = qw2qt(s, p.lambda);
  set(St, 'Error',   0);
  set(St, 'Monitor', 1); 
  St   = St.*q./Kf;   % Schober Eq (9.294)
  
  clear Ef Vf Kf
  % resample histogram when axes are not vectors
  hist_me = false;
  for index=1:ndims(St)
    x = getaxis(St, index);
    if numel(x) ~= length(x), hist_me=true; break; end
  end
  if hist_me
    St   = hist(St, size(St));  % the data set is well covered here and interpolation works well.
  end
  St   = trapz(St, 2);
  
  % Elastic peak position
  EPP = Sqw_getT(s, {'ElasticPeakPosition' 'Elastic_peak_channel' 'Elastic_peak' 'Peak_channel' 'Elastic'});
  if isempty(EPP)
    St  = rmaxis(St, 1); % to channels
    EPP = Sqw_getEPP(St);
  end
  
  St = double(St);
  
  % reduce St to max 512 elements (limit in MUPHOCOR)
  if numel(St) > 512
    x  = 1:numel(St);
    nx = linspace(1,numel(St),512);
    EPP= EPP*512/numel(St);
    St = interp1(x, St, nx);
  end

  % we get the angle range
  phi = getaxis(qw2phiw(s, p.lambda), 2);
  
  % This is equivalent to: S(q,w) -> S(phi,channel)
  %   [~,spc] = qw2phi(s, p.lambda);           
  %   St   = trapz(spc .* sind(phi),2);  is very slow in interpolation
  % This computation corresponds with (when using a LAMP/NeXus data set). 'Y' is angle; 
  %   St = sum(f.DATA.*repmat(sind(f.Y)',[size(f.DATA,1),1]),2);
  % except for a normalisation factor

  disp([ mfilename ': Computing vDOS for ' s.Tag ' ' s.Title ' from ' s.Source ]);
  if isfield(s, 'ChemicalFormula')
    disp([ '  Material     ' get(s, 'ChemicalFormula') ]);
  else
    disp([ '  Material     ' char(s) ]);
  end
  disp([ '  Temperature= ' num2str(p.temp0)  ' [K]' ])
  disp([ '  Wavelength = ' num2str(p.lambda)  ' [Angs] incident' ])
  disp([ '  Q          = [' num2str(min(q(:))) ' ' num2str(max(q(:))) '] [Angs-1]' ])
  disp([ '  w          = [' num2str(min(w(:))) ' ' num2str(max(w(:))) '] [meV]' ])
  disp([ '  theta      = [' num2str(min(phi(:))) ' ' num2str(max(phi(:))) '] [deg]' ])
  
  % transfer parameters from iData_Sqw2D to muphocor input file (as a structure)
  % nspec can be > 1
  % then: amasi, conci and sigi must be vectors of length nspec
  params = muphocor_create_params_struc( ...
    'fp',     p.distance*100.0, ...
    'fimi',   min(phi(:)), ...
    'fima',   max(phi(:)), ...
    'hox',    max(abs(w(:))), ...
    'nspec',  numel(p.sigi), ...
    'no',     numel(St), ...
    'noo',    numel(St), ...
    'nu',     1, ...
    'unt',    min(St(:)), ...
    'nuu',    1);
  
  clear q w
  
  params.xnel  = EPP;
  params.lm    = min(numel(St), 512);  % LM must be limited to 512 (LM2 and LM5 to 1024 in MUPHOCOR)
  params.title = char(s);
  
  % transfer any additional parameter from input arguments to MUPHOCOR
  for f=fieldnames(p)'
    params.(f{1}) = p.(f{1});
  end

  % run MUPHOCOR
  [g, params.output, params.input] = muphocor_run(compiled, params, St);
  
  % g is an array with columns:
  % [index] [hw] [ g(w) g_0(w) g_multi(w) ]
  hw = g(:,2);
  
  gw = iData(hw, g(:,3)); 
  gw.Label ='gDOS';
  gw.Title = [ gw.Label '(' title(s) ')' ];
  title(gw, [ gw.Label ' g(w) [meV^{-1}]' ]); 
  
  g0w= iData(hw, g(:,4));
  g0w.Label ='vDOS';
  g0w.Title = [ g0w.Label '(' title(s) ')' ];
  title(g0w, [ g0w.Label ' g(w) [meV^{-1}]' ]); 
  
  gmw= iData(hw, g(:,5));
  gmw.Label ='gDOS multi-phonons';
  gmw.Title = [ gmw.Label '(' title(s) ')' ];
  title(gmw, [ gmw.Label ' g(w) [meV^{-1}]' ]); 
  
  g = [ gw g0w gmw ];
  xlabel(g, 'Energy [meV]'); set(g, 'Error', 0, 'Monitor', 1);
  set(g, 'UserData', s.UserData);
  setalias(g, 'muphocor',   params);
  
  f = getalias(s);
  for index=1:numel(getalias(s))
    if ~isfield(g, f{index})
      [link, lab] = getalias(s, f{index});
      g = setalias(g, f{index}, link, lab);
    end
  end

end % muphocor

% ------------------------------------------------------------------------------
function muphocor_write_input(foutname, par, St)
  % Write the input.txt file
  % input parameters, then the sin(theta) S(theta,t)
  [fid, mess] = fopen(foutname, 'w');
  if fid == -1
    disp(mess)
    error([ mfilename ': could not create the MUPHOCOR input file ' foutname ]); 
  end
  
  fprintf(fid, '%s\n', par.title);
  fmt = [repmat('%10.2f',1,8), '\n'];
  fprintf(fid, fmt, par.e0, par.fp, par.cw, par.xnel, par.fimi, par.fima, par.abk, par.aemp);
  fmt = [repmat('%5d',1,4), '\n'];
  fprintf(fid, fmt, par.lm, par.nspec, par.iquo, par.lpar);
  fmt = [repmat('%10.2f',1,3), '\n'];
  fprintf(fid, fmt, par.hox, par.temp0, par.dw);
  for k=1:par.nspec
      fprintf(fid, fmt, par.amasi(k), par.sigi(k), par.conci(k));
  end
  fmt = [repmat('%5d',1,8), '\n'];
  fprintf(fid, fmt, par.npho, par.itm, par.ivit, par.idw, par.iemp, par.ires, par.ipr, par.iloss);
  fmt = [repmat('%5d',1,6), '\n'];
  fprintf(fid, fmt, par.isour, par.nuu, par.noo, par.nu, par.no, par.iglu);
  fmt = [repmat('%10.4f',1,2), '\n'];
  fprintf(fid, fmt, par.unt, par.fun);
  fmt = [repmat('%10.4f',1,8), '\n'];
  fprintf(fid, fmt, St);
  fmt = '\n%10.4f';
  fprintf(fid, fmt, 0.0);

  fclose(fid);
end % muphocor_write_input

function [s, mess, in] = muphocor_run(compiled, params, St)
  % run MUPHOCOR for given parameters and channel spectrum [\int sin(theta) S(theta, channel) dtheta]
  
  % create a temporary directory to work in
  tmp = tempname;
  mkdir(tmp);
  mess = []; s = [];
  
  % write input muphocor file
  muphocor_write_input(fullfile(tmp, 'input.txt'), params, St);
  in = fileread(fullfile(tmp, 'input.txt'));
  
  % run MUPHOCOR executable
  if ismac,      precmd = 'DYLD_LIBRARY_PATH= ; DISPLAY= ; ';
  else0f isunix, precmd = 'LD_LIBRARY_PATH= ;  DISPLAY= ; '; 
  else           precmd = ''; end

  % If defined , go directly to:
  cmd = [ compiled ' < ' fullfile(tmp, 'input.txt') ];

  % Perform the calculation
  disp(cmd);
  p = pwd;
  cd(tmp);  % go in tmp dir
  try
    [status, mess] = system([ precmd cmd ]);
  end
  
  % read results
  cd(p);  % revert to pwd
  
  if exist(fullfile(tmp, 'GDOS_plot_f90.dat'), 'file')
    s = load(fullfile(tmp, 'GDOS_plot_f90.dat'));
  else
    disp(mess);
    error([ mfilename ': ERROR: Calculation failed. Check parameters and output above.' ])
  end
  
end
