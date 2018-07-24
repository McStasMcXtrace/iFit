function g = muphocor(s, lambda, T, mass, sigma_inc, conc)
  % run MUPHOCOR on the given S(q,w) data set, incident wavelength and Temperature
  %
  %   g = muphocor(s, lambda, T, mass, sigma_inc, conc)
  %
  % When not given, the physical parameters are searched in the object.
  % The material physical properties include the masses, the incoherent 
  %   scattering cross sections, and the concentrations which can be given as 
  %   vectors per element/non equivalent atom in the material.
  % The concentrations are the relative proportion of atoms in the material.
  %
  % input:
  %   s:      iData_Sqw2D object S(q,w)
  %   lambda: optional incident neutron wavelength [Angs]
  %   T:      optional Temperature [K]
  %   mass:   optional material mass [g/mol]
  %   sigma_inc: optional incoherent scattering cross section [barns]
  %   conc:   optional material concentration [1]
  
  % check if muphocor is compiled, else compile it
  persistent compiled
  
  if isempty(compiled)
    compiled = muphocor_compile_binary; % check and possibly compile MUPHOCOR. Return path the EXE
  end
  if nargin < 2, lambda   =[]; end
  if nargin < 3, T        =[]; end
  if nargin < 4, mass     =[]; end
  if nargin < 5, sigma_inc=[]; end
  if nargin < 6, conc     =[]; end
  
  if isempty(lambda)
    [s,lambda,distance,chwidth,~,wavevector] = Sqw_search_lambda(s);
  else
    [s,~,distance,chwidth,~,wavevector]      = Sqw_search_lambda(s);
  end
  Ei    = 81.805./lambda.^2;
  if chwidth < 1e-4
    chwidth = chwidth*1e6;  % [s] to [us];
  end
  
  if isempty(T)
    T    = Sqw_getT(s);
  end
  if isempty(T)
    T = 293;
    disp([ mfilename ': WARNING: Using Temperature=' num2str(T) ' [K] for data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
  end
  
  % check for other required parameters: mass/weight, b_inc, concentration
  if isempty(mass)
    mass = Sqw_getT(s, {'weight'	'mass' 'AWR'});
  end
  if isempty(mass)
    disp( [ mfilename ': Set object ' inputname(1) '.mass to the total molar weight, or a vector with molar weight per atom in the material.' ]);
    error([ mfilename ': ERROR: Unspecified mass [g/mol] in data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
  end
  if isempty(sigma_inc)
    sigma_inc = Sqw_getT(s, {'sigma_inc'});
  end
  if isempty(sigma_inc)
    disp( [ mfilename ': Set object ' inputname(1) '.sigma_inc to the total incoherent neutron scattering cross section, or a vector with sigma_inc per atom in the material.' ]);
    error([ mfilename ': ERROR: Unspecified sigma_inc [barn] in data set ' s.Tag ' ' s.Title ' from ' s.Source ]);
  end
  if isempty(conc)
    conc = Sqw_getT(s, {'concentration'});
  end
  if isempty(conc)
    conc = ones(size(mass));
  end
  
  % data must be [\int sin(theta) S(theta,t) dtheta] that is start from time distribution
  [~,spc] = qw2phi(s, lambda);  % S(q,w) -> S(phi,channel)
  phi  = getaxis(spc, 2);
  w    = getaxis(s, 1);
  q    = getaxis(s, 2);

  
  
  disp([ mfilename ': Computing vDOS for ' s.Tag ' ' s.Title ' from ' s.Source ]);
  disp([ '  Temperature= ' num2str(T)  ' [K]' ])
  disp([ '  Wavelength = ' num2str(lambda)  ' [Angs]' ])
  disp([ '  Q          = [' num2str(min(q(:))) ' ' num2str(max(q(:))) '] [Angs-1]' ])
  disp([ '  w          = [' num2str(min(w(:))) ' ' num2str(max(w(:))) '] [meV]' ])
  disp([ '  theta      = [' num2str(min(phi(:))) ' ' num2str(max(phi(:))) '] [deg]' ])
  
  St   = trapz(spc .* sind(phi), 2);
  St   = double(St);
  
  % this computation corresponds with (when using a LAMP/NeXus data set). 'Y' is angle; 
  %   St = sum(f.DATA.*repmat(sind(f.Y)',[size(f.DATA,1),1]),2);
  % except for a normalisation factor

  % transfer parameters from iData_Sqw2D to muphocor input file (as a structure)
  % nspec can be > 1
  % then: amasi, conci and sigi must be vectors of length nspec
  params = muphocor_create_params_struc( ...
    'e0',     Ei, ...
    'fp',     distance*100.0, ...
    'cw',     chwidth, ...
    'fimi',   min(phi(:)), ...
    'fima',   max(phi(:)), ...
    'temp0',  T, ...
    'hox',    max(w(:)), ...
    'nspec',  1, ...
    'amasi',  mass,  ...
    'conci',  1, ...
    'sigi',   sigma, ...
    'no',     numel(St), ...
    'noo',    numel(St), ...
    'nu',     1, ...
    'nuu',    1);
  params.title = char(s);
  
  % Elastic peak position
  params.xnel = Sqw_getEPP(spc);
  params.lm   = min(numel(St), 512);  % LM must be limited to 512 (LM2 and LM5 to 1024 in MUPHOCOR)
  
  % run MUPHOCOR
  g = muphocor_run(compiled, parameters, St);
  
  % g is an array with columns:
  % [index] [hw] [ g(w) g_0(w) g_multi(w) ]

end % muphocor

% ------------------------------------------------------------------------------
function muphocor_write_input(foutname, params, St)
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

function s = muphocor_run(compiled, params, St)
  % run MUPHOCOR for given parameters and channel spectrum [\int sin(theta) S(theta, channel) dtheta]
  
  % create a temporary directory to work in
  tmp = tempname;
  mkdir(tmp);
  
  % write input muphocor file
  muphocor_write_input(fullfile(tmp, 'input.txt'), params, St);
  
  % run MUPHOCOR executable
  if ismac,      precmd = 'DYLD_LIBRARY_PATH= ; DISPLAY= ; ';
  elseif isunix, precmd = 'LD_LIBRARY_PATH= ;  DISPLAY= ; '; 
  else           precmd = ''; end

  % If defined , go directly to:
  cmd = [ compiled ' < ' fullfile(tmp, 'input.txt') ];

  % Perform the calculation
  disp(cmd);
  p = pwd;
  cd(tmp);  % go in tmp dir
  try
    [status, mess] = system([ precmd cmd ]);
    disp(mess)
  end
  
  % read results
  cd(p);  % revert to pwd
  
  s = load(fullfile(tmp, 'GDOS_plot_f90.dat'));
  
end
