function options=sqw_phonons_argin(varargin)
% sqw_phonons_argin: extracts options from the arguments
%
% returns an 'options' structure.

% defaults
options.supercell  = 4;
options.calculator = '';
options.kpoints    = 4;
options.xc         = 'PBE';
options.mode       = 'fd';            % GPAW
options.potentials = '';
options.diagonalization = ''; % GPAW would prefer rmm-diis
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
    elseif strcmpi(varargin{index},'qe') || strcmpi(varargin{index},'espresso') || strcmpi(varargin{index},'quantumespresso')
      options.calculator = 'quantumespresso';
    elseif strcmpi(varargin{index},'autoplot')
      options.autoplot = 1;
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
if ~isfield(options,'target')
  options.target = tempname; % everything will go there
  mkdir(options.target)
end
if isfield(options, 'smearing') && isempty(options.occupations)
  options.occupations = options.smearing;
end
if isscalar(options.supercell)
  options.supercell=[ options.supercell options.supercell options.supercell ]; 
end
if isscalar(options.kpoints)
  options.kpoints=[ options.kpoints options.kpoints options.kpoints ]; 
end
