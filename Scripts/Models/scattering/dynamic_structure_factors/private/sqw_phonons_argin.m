function options=sqw_phonons_argin(varargin)
% sqw_phonons_argin: extracts options from the arguments
%
% returns an 'options' structure.

% defaults
options.supercell  = 3;
options.calculator = '';
options.kpoints    = 3;
options.xc         = 'PBE';
options.mode       = 'fd';            % GPAW
options.potentials = '';
options.diagonalization = 'rmm-diis'; % GPAW would prefer rmm-diis
options.occupations= '';
options.ecut       = 0;
options.nbands     = 0;
options.nsteps     = 0;
options.toldfe     = 0;
options.command    = '';
options.raw        = '';
options.autoplot   = 0;
options.gui        = 0;
options.htmlreport = 0;
options.dos        = 0;
options.email      = '';

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
    elseif strcmpi(varargin{index},'semiconductor')
      options.occupations = 'semiconductor';
    elseif strcmpi(varargin{index},'dos')
      options.dos = 1;
    elseif strcmpi(varargin{index},'plot') || strcmpi(varargin{index},'autoplot')
      options.autoplot = 1;
    elseif strcmpi(varargin{index},'html') || strcmpi(varargin{index},'htmlreport') || strcmpi(varargin{index},'report')
      options.htmlreport = 1;
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
    elseif strcmpi(varargin{index},'qe') || strcmpi(varargin{index},'espresso') || strcmpi(varargin{index},'quantumespresso')
      options.calculator = 'quantumespresso';
    elseif strcmpi(varargin{index},'gui')
      options.gui = 1;
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
if isfield(options,'dir') && ~isfield(options,'target')
  options.target = options.dir;
end
if ~isfield(options,'target')
  options.target = tempname; % everything will go there
end
if isfield(options, 'smearing') && isempty(options.occupations)
  options.occupations = options.smearing;
end
if isfield(options, 'kpts')
  options.kpoints=options.kpts;
end
if isfield(options, 'plot')
  options.autoplot=options.plot;
end
if isfield(options, 'html')
  options.htmlreport=options.html; options.dos=1;
end
if isscalar(options.supercell)
  options.supercell=[ options.supercell options.supercell options.supercell ]; 
end
if isscalar(options.kpoints)
  options.kpoints=[ options.kpoints options.kpoints options.kpoints ]; 
end
if options.htmlreport == 1, options.dos=1; end

% make sure target is a fully qualified path
if options.target(1) ~= filesep
  options.target = fullfile(pwd, options.target);
end
if ~isdir(options.target)
  mkdir(options.target);
end

