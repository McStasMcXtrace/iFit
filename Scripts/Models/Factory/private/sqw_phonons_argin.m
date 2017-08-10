function options=sqw_phonons_argin(varargin)
% sqw_phonons_argin: extracts options from the arguments
%
% returns an 'options' structure.

% defaults
options.supercell  = 0; % automatic
options.calculator = '';
options.kpoints    = 0; % automatic
options.xc         = 'PBE';
options.mode       = 'pw';            % GPAW
options.potentials = '';
options.diagonalization = 'rmm-diis'; % GPAW would prefer rmm-diis, also fastest
options.occupations= 'auto';          % could be auto, fixed, metal, semiconductor
options.ecut       = 0;
options.nbands     = 0;
options.nsteps     = 0;
options.toldfe     = 1e-5;  % required. When too low, the convergence gets difficult.
options.command    = '';
options.raw        = '';
options.autoplot   = 0;
options.gui        = nan;
options.htmlreport = 0;
options.optimizer  = '';
options.accuracy   = 'very fast';  % can be 'fast' 'very fast' or 'accurate' (much slower)
options.disp       = 0.01;    % displacement of atoms in Angs
options.use_phonopy= 1;

% read input arguments
for index=1:numel(varargin)
  if ischar(varargin{index}) && isempty(dir(varargin{index})) && ~isempty(find(varargin{index} == '='))
    % first try to build a structure from the string
    this = str2struct(varargin{index});
    if isstruct(this)
      varargin{index} = this;
    end
  end
  if ischar(varargin{index})
    [p,f,e] = fileparts(varargin{index});
    % handle static options: metal,insulator, random
    if strcmpi(varargin{index},'smearing') || strcmpi(varargin{index},'metal')
      options.occupations = 'smearing';
    elseif strcmpi(varargin{index},'fixed') || strcmpi(varargin{index},'insulator')
      options.occupations = 'fixed';
    elseif strcmpi(varargin{index},'semiconductor') || strcmpi(varargin{index},'semi')
      options.occupations = 'semiconductor';
    
      
% calculators
    elseif strcmpi(varargin{index},'emt')
      options.calculator = 'EMT';
    elseif strcmpi(varargin{index},'gpaw')
      options.calculator = 'GPAW';
    elseif strcmpi(varargin{index},'jacapo') || strcmpi(varargin{index},'dacapo')
      options.calculator = 'Jacapo';
    elseif strcmpi(varargin{index},'nwchem')
      options.calculator = 'NWChem';
    elseif strcmpi(varargin{index},'elk')
      options.calculator = 'Elk';
    elseif strcmpi(varargin{index},'abinit')
      options.calculator = 'ABINIT';
    elseif strcmpi(varargin{index},'vasp')
      options.calculator = 'VASP';
    elseif strcmpi(varargin{index},'qe_ase') || strcmpi(varargin{index},'espresso_ase') ...
      || strcmpi(varargin{index},'quantumespresso_ase')
      options.calculator = 'quantumespresso_ase';
    elseif strcmpi(varargin{index},'qe_phon') || strcmpi(varargin{index},'espresso_phon') ...
      || strcmpi(varargin{index},'quantumespresso_phon') || strcmpi(varargin{index},'phon') ...
      || strcmpi(varargin{index},'qe') || strcmpi(varargin{index},'espresso') ...
      || strcmpi(varargin{index},'quantumespresso')
      options.calculator = 'quantumespresso';
    elseif strcmpi(varargin{index},'octopus')
      options.calculator = 'octopus';
    
    % other options
    elseif strcmpi(varargin{index},'plot') || strcmpi(varargin{index},'autoplot')
      options.autoplot = 1;
    elseif strcmpi(varargin{index},'html') || strcmpi(varargin{index},'htmlreport') || strcmpi(varargin{index},'report')
      options.htmlreport = 1;
    elseif strcmpi(varargin{index},'gui')
      options.gui = 'init';
    elseif strcmpi(varargin{index},'optimize') || strcmpi(varargin{index},'optimise') || strcmpi(varargin{index},'minimize')
      options.optimizer = 'BFGS';
    elseif strcmpi(varargin{index},'fast') || strcmpi(varargin{index},'low') || strcmpi(varargin{index},'coarse')
      options.accuracy = 'fast';
      options.use_phonopy = 0;
    elseif strcmpi(varargin{index},'very fast')
      options.accuracy = 'very fast';
    elseif strcmpi(varargin{index},'accurate') || strcmpi(varargin{index},'high') || strcmpi(varargin{index},'slow')
      options.accuracy = 'accurate';
      options.use_phonopy = 0;
    elseif strcmpi(varargin{index},'phonopy') 
      options.use_phonopy = 1;
      options.accuracy = 'very fast';
    end
  end
  if isstruct(varargin{index})
    % a structure: we copy the fields into options.
    this = varargin{index};
    f    =fieldnames(this);
    for i=1:numel(fieldnames(this))
      options.(f{i}) = this.(f{i});
    end
  end
end
if isfield(options, 'optimize')
  options.optimizer=options.optimize;
end
if isfield(options,'dir') && ~isfield(options,'target')
  options.target = options.dir;
end
if ~isfield(options,'target')
  options.target = tempname; % everything will go there
end
if isfield(options, 'smearing') && isempty(options.occupations)
  options.occupations = options.smearing;
end
if isfield(options, 'nstep') && ~options.nsteps
  options.nsteps = options.nstep;
end
if isfield(options, 'cutoff') && ~options.ecut
  options.ecut = options.cutoff;
end
if isfield(options, 'kpts')
  options.kpoints=options.kpts;
end
if isfield(options, 'plot')
  options.autoplot=options.plot;
end
if isfield(options, 'html')
  options.htmlreport=options.html;
end
if isscalar(options.supercell)
  options.supercell=[ options.supercell options.supercell options.supercell ]; 
end
if isscalar(options.kpoints)
  options.kpoints=[ options.kpoints options.kpoints options.kpoints ]; 
end

% make sure target is a fully qualified path
if options.target(1) ~= filesep
  options.target = fullfile(pwd, options.target);
end
if ~isdir(options.target)
  mkdir(options.target);
end

