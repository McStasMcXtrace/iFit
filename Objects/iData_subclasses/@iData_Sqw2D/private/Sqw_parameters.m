function  [s,parameters,fields] = Sqw_parameters(s, fields)
% [s,p,fields] = Sqw_parameters(s,type): search for parameter values in a S(q,w) data set
%
% iData_Sqw2D: Sqw_parameters: search for physical quantities in object.
%   This search is also done when creating iData_Sqw2D objects.
%
% input:
%   s: any iData object, including S(q,w) and S(alpha,beta) ones.
%   fields: an optional list of items to search (cellstr)
% output:
%   s: updated object with found parameters
%   p: parameters as a structure, also stored into s.parameters
%
% (c) E.Farhi, ILL. License: EUPL.

% extract a list of parameters
parameters = [];

if nargin == 0, return; end
if ~isa(s, 'iData'), s=iData(s); end
if nargin < 2, fields=[]; end
if isempty(fields) % default: search all parameters and assemble them
  [s,parameters1,fields1] = Sqw_parameters(s, 'sqw');
  [s,parameters2,fields2] = Sqw_parameters(s, 'sab');
  if isstruct(parameters1) && isstruct(parameters2)
    % merge the structures/cells handling duplicated fields
    f = [ fieldnames(parameters1) ; fieldnames(parameters2) ];
    [pnames,index] = unique(f);
    pairs = [fieldnames(parameters1), struct2cell(parameters1); ...
             fieldnames(parameters2), struct2cell(parameters2)].';
    parameters = cell2struct(pairs(2,index), pairs(1,index), 2);
    fields = [ fields1 fields2 ];
    fields = fields(index);
    setalias(s, 'parameters', parameters(1));
  end
  if nargout == 0 & length(inputname(1))
    assignin('caller',inputname(1),s);
  end
  return
end

if numel(s) > 1
  for index=1:numel(s)
    [s(index),parameters,fields] = Sqw_parameters(s(index), fields);
  end
  return
end

% make sure we have a 'parameters' field
if isfield(s, 'parameters')
  parameters = get(s, 'parameters');
end

if ischar(fields) && strcmpi(fields, 'sqw')
  % a list of parameter names, followed by comments
  % aliases can be specified when a parameter is given as a cellstr, with 
  % consecutive possibilities.
  fields={ ...
      'density [g/cm3] Material density'	...
     {'weight [g/mol] Material molar weight'	'molar_mass' 'mass' 'masses' 'AWR'}...
      'T_m [K] Melting T'	...
      'T_b [K] Boiling T'	'MD_at Number of atoms in the molecular dynamics simulation'	...
      'MD_box [Angs] Simulation box size'	...
     {'Temperature [K]' 'T' 'TEMP'}	...
      'dT [K] T accuracy' ...
      'MD_duration [ps] Molecular dynamics duration'	'D [cm^2/s] Diffusion coefficient' ...
      'sigma_coh [barn] Coherent scattering neutron cross section' ...
      'sigma_inc [barn] Incoherent scattering neutron cross section'	...
      'b_coh [fm] Coherent scattering length (neutron)' ...
      'b_inc [fm] Incoherent scattering length (neutron)'	...
      'sigma_abs [barn] Absorption neutron cross section'	...
      'c_sound [m/s] Sound velocity' ...
      'At_number Atomic number Z' ...
      'Pressure [bar] Material pressure' ...
      'v_sound [m/s] Sound velocity' ...
      'Material Material description' ...
      'Phase Material state' ...
      'Scattering process' ...
     {'Wavelength [Angs] neutron wavelength' 'lambda' 'Lambda' } ...
     {'Instrument neutron spectrometer used for measurement' 'SOURCE' } ...
      'Comment' ...
      'C_s [J/mol/K]   heat capacity' ...
      'Eta [Pa s] viscosity' ...
      'L [J/mol] latent heat' ...
     {'classical [0=from measurement, with Bose factor included, 1=from MD, symmetric]' 'symmetry' } ...
      'multiplicity  [atoms/unit cell] number of atoms/molecules per scattering unit cell' ...
     {'IncidentEnergy [meV] neutron incident energy' 'fixed_energy' 'energy' 'ei' 'Ei'} ...
     {'FinalEnergy [meV] neutron final energy' 'ef' 'Ef'} ...
     {'IncidentWavevector [Angs-1] neutron incident wavevector' 'ki' 'Ki'} ...
     {'FinalWavevector [Angs-1] neutron final wavevector' 'kf' 'Kf'} ...
     {'Distance [m] Sample-Detector distance' 'distance' } ...
     {'ChannelWidth [time unit] ToF Channel Width' 'Channel_Width'} ...
     {'V_rho [Angs^-3] atom density' 'rho' } ...
     {'DetectionAngles [deg] Detection Angles'} ...
     {'ElasticPeakPosition [chan] Channel for the elastic peak' 'Elastic_peak_channel' 'Elastic_peak' 'Peak_channel' 'Elastic'} ...
     {'ChemicalFormula Material Chemical Formula' 'chemical_formula' 'formula' 'chemical' 'ZSYMAM'} ...
     {'Concentration Proportional concentration of atoms in material'} ...
    };
elseif ischar(fields) && strcmpi(fields, 'sab')
  % ENDF compatibility flags: MF1 MT451 and MF7 MT2/4
  % we have added parameters from PyNE/ENDF
  fields = { ...
    'MF ENDF File number' ...
    'MT ENDF Reaction type' ...
    'ZA ENDF material id (Z,A)' ...
   {'weight [g/mol] Material molar weight'	'mass' 'AWR'} ...
   {'classical [0=from measurement, with Bose factor included, 1=from MD, symmetric]' 'symmetry' } ...
   {'Temperature [K]' 'T' }	...
    'sigma_coh [barn] Coherent scattering neutron cross section' ...
    'sigma_inc [barn] Incoherent scattering neutron cross section'	...
    'sigma_abs [barn] Absorption neutron cross section'	...
    'b_coh [fm] Coherent scattering length (neutron)' ...
    'b_inc [fm] Incoherent scattering length (neutron)'	...
    'At_number Atomic number Z' ...
    'Scattering process' ...
   {'Material ENDF material code' 'MAT' } ...
    'charge' ...
    'EDATE ENDF Evaluation date' ...
    'LRP ENDF Resonance flag' ...
   {'LFI ENDF Fission flag' 'fissionable' } ...
   {'NLIB ENDF Library' 'library' } ...
   {'NMOD ENDF modification id' 'modification' } ...
   {'ELIS ENDF Excitation energy' 'excitation_energy' } ...
   {'STA ENDF Target stability flag' 'stable' } ...
   {'LIS ENDF State number of the target nucleus' 'state' } ...
   {'LISO ENDF Isomeric state number' 'isomeric_state' } ...
    'NFOR ENDF format id' ...
    'AWI ENDF Projectile mass in neutron units' ...
   {'EMAX ENDF Max energy [eV]' 'energy_max' } ...
    'LREL ENDF Release number' ...
   {'NSUB ENDF Sublibrary' 'sublibrary' } ...
    'NVER ENDF Release number' ...
    'TEMP ENDF Target T for Doppler broadening [K]' ...
   {'LDRV ENDF derived evaluation id' 'derived' }...
    'NWD ENDF Number of records' ...
    'NXC ENDF Number of records in directory' ...
   {'ZSYMAM ENDF Character representation of the material' 'ChemicalFormula' 'chemical_formula' 'formula'} ...
   {'ALAB ENDF Laboratory' 'laboratory' } ...
   {'AUTH ENDF Author' 'author' } ...
   {'REF ENDF Primary reference' 'reference' } ...
   {'DDATE ENDF Distribution date' 'date_distribution' } ...
   {'RDATE ENDF Release date' 'date_release' } ...
   {'ENDATE ENDF Entry date' 'date_entry' } ...
    'HSUB ENDF Libabry hierarchy' ...
    'SB ENDF bound cross section [barns]' ...
    'NR' ...
   {'NP ENDF Nb of pairs' 'n_pairs' } ...
    'LTHR ENDF Thermal data flag' ...
    'S ENDF Thermal elastic' 'E' 'LT ENDF T dependance flag' ...
   {'LAT ENDF T used for (alpha,beta)' 'temperature_used' } ...
    'LASYM ENDF Asymmetry flag' 'B ENDF model parameters' 'NI' ...
   {'NS ENDF Number of states/parameters' 'num_non_principal' } ...
   {'LLN ENDF storage form [linear/ln(S)]' 'ln_S_' } ...
    'NT' 'Teff ENDF Effective T [K]' ...
   {'Sab S(alpha,beta,T) scattering law' 'scattering_law' }
        };
end

% add parameters by searching field names from 'fields' in the object
for index=1:length(fields)
  f = fields{index};
  if ischar(f), f= {f}; end
  p_name = strtok(f{1});          % the first name before the comment. We will assign it.
  for f_index=1:numel(f)
    name   = strtok(f{f_index});  % the name of the field or alternatives
    if     isfield(s, p_name) val0=get(s, p_name);
    elseif isfield(parameters, name)
      val0 = parameters.(name);
    elseif isfield(parameters, p_name)
      val0 = parameters.(p_name);  
    elseif isfield(s, name)   val0=get(s, name);
    else val0=[]; end
    if index==1 && f_index==1
      links    = findfield(s,strtok(f{f_index}),'exact case');
      if isempty(links)
        links    = findfield(s,strtok(f{f_index}),'exact');
      end
    else
      links    = findfield(s,strtok(f{f_index}),'exact cache case');
      if isempty(links)
        links    = findfield(s,strtok(f{f_index}),'exact cache');
      end
      if isempty(links) && numel(strtok(f{f_index})) > 3
        links    = findfield(s,strtok(f{f_index}),'cache'); % search for incomplete names
      end
    end
    if isfield(s.Data, name)
      links = [ 'Data.' name ];
    elseif isfield(parameters, name)
      links = ['parameters.' name ];
    end
    if ~iscell(links), links = { links }; end
    for l_index=1:numel(links)
      link = links{l_index};
      try
        val = get(s,link);
      catch
        link = []; val=[];
      end
      % Do not update when:
      %   old value is 'larger' than previous one
      if ~isempty(val) && ~isempty(val0)
        if (isnumeric(val) && isnumeric(val0) && isscalar(val) && isscalar(val0) && val < val0) ...
        || (isnumeric(val) && isnumeric(val0) && numel(val) < numel(val0))
          link=[]; 
        end
      end
      % make sure it makes sense to update the link. must not link to itself.
      if ~isempty(link) && ~isempty(val) && ~isfield(parameters, p_name)
        if numel(val) < 100 || strcmp(p_name, link)
          parameters.(p_name) = val;   % update/store the value
        else 
          parameters.(p_name) = link;  % update/store the link/alias
        end
      end
    end

  end
  fields{index} = f{1}; % store the parameter name, plus optional label
end

% check for chemical formula
[s, parameters] = Sqw_parameters_chemicalformula(s, parameters);

% update parameter values (replace links by values) and object aliases
for index=1:numel(fields)
  [name, comment] = strtok(fields{index});
  if isfield(parameters, name)
    link = parameters.(name);  % link
    % get and store the value of the parameter
    if ~isempty(link) && ~isempty(val) && (~ischar(link) || ~strcmp(name, link))  % must not link to itself
      s=setalias(s, name, [ 'parameters.' name ], strtrim(comment));  % link to parameters
      try
        parameters.(name) = get(s, link); % the value
      catch
        if isfield(s, 'name')
          parameters.(name) = get(s, name); % the value
        end
      end
    end
  end
end

% now transfer parameters into the object, as alias values
s=setalias(s, 'parameters',       parameters, 'Material parameters');

if nargout == 0 & length(inputname(1))
  assignin('caller',inputname(1),s);
end


% ------------------------------------------------------------------------------
function [s, parameters] = Sqw_parameters_chemicalformula(s, parameters)
  % Sqw_parameters_chemicalformula: search for a chemical formula, and derive 
  % missing mass, scattering cross sections
  
  Z = {'H' 'He' 'Li' 'Be' 'B' 'C' 'N' 'O' 'F' 'Ne' 'Na' 'Mg' 'Al' 'Si' 'P' 'S' 'Cl' 'Ar' 'K' 'Ca' 'Sc' 'Ti' 'V' 'Cr' 'Mn' 'Fe' 'Co' 'Ni' 'Cu' 'Zn' 'Ga' 'Ge' 'As' 'Se' 'Br' 'Kr' 'Rb' 'Sr' 'Y' 'Zr' 'Nb' 'Mo' 'Tc' 'Ru' 'Rh' 'Pd' 'Ag' 'Cd' 'In' 'Sn' 'Sb' 'Te' 'I' 'Xe' 'Cs' 'Ba' 'La' 'Ce' 'Pr' 'Nd' 'Pm' 'Sm' 'Eu' 'Gd' 'Tb' 'Dy' 'Ho' 'Er' 'Tm' 'Yb' 'Lu' 'Hf' 'Ta' 'W' 'Re' 'Os' 'Ir' 'Pt' 'Au' 'Hg' 'Tl' 'Pb' 'Bi' 'Po' 'At' 'Rn' 'Fr' 'Ra' 'Ac' 'Th' 'Pa' 'U' 'Np' 'Pu' 'Am' 'Cm' 'Bk' 'Cf' 'Es' 'Fm' 'Md' 'No' 'Lr' 'Rf' 'Db' 'Sg' 'Bh' 'Hs' 'Mt' 'Ds' 'Rg' 'Cn' 'Nh' 'Fl' 'Mc' 'Lv' 'Ts' 'Og'};
  
  if isfield(s, 'ChemicalFormula') && ischar(get(s,'ChemicalFormula')) ...
    && ~isempty(get(s,'ChemicalFormula')) ...
    && ~strcmp(getalias(s, 'ChemicalFormula'), 'parameters.ChemicalFormula')
    parameters.ChemicalFormula = get(s,'ChemicalFormula');
  end
  
  if (isfield(parameters, 'At_number') && isnumeric(parameters.At_number)) ...
    && (~isfield(s, 'ChemicalFormula') || ~ischar(get(s,'ChemicalFormula'))) ...
    && (~isfield(parameters, 'ChemicalFormula') || ~ischar(parameters.ChemicalFormula)) ...
    && all(parameters.At_number > 0)
    parameters.ChemicalFormula = sprintf('%s ', Z{parameters.At_number});
  end

  if ~isfield(parameters, 'ChemicalFormula') || ~ischar(parameters.ChemicalFormula) 
    at = 'A[cglmrstu]|B[aehikr]?|C[adeflmnorsu]?|D[bsy]|E[rsu]|F[elmr]?|G[ade]|H[efgos]?|I[nr]?|Kr?|L[airuv]|M[dgnot]|N[abdeiop]?|Os?|P[abdmortu]?|R[abefghnu]|S[bcegimnr]?|T[abcehilm]|U(u[opst])?|V|W|Xe|Yb?|Z[nr]';
    qt ='\d*\.?\d*';
    % regexp to match formula with only Atom symbols
    % r='[A-Z][a-z]?\d*\.?\d*|(?<!\([^)]*)\(.*\)\d+\.?\d*(?![^(]*\))'      also works
    % r='\(?([A-Z]{1}[a-z]?[^a-z])(\d*\.?\d*)([+-]?)(\d*\.?\d*)\)?(\d*)';  works nicely, but not atom specific
    r= [ '\(?(' at ')([^a-z\._-])(' qt ')([+-]?)(' qt ')\)?(\d*)' ];
    form = regexp(findstr(s), r, 'match');
    % get the longest formula
    try
      [~,index] = max(cellfun(@numel, form));
      form = form{index};
      parameters.ChemicalFormula = sprintf('%s', form{:});
    end
  end

  % determine mass, and neutron cross sections
  % we parse the chemical formula and split it into tokens
  if isfield(parameters, 'ChemicalFormula') && ischar(parameters.ChemicalFormula)
    try
      r = parse_formula(parameters.ChemicalFormula);
      if ~isfield(parameters, 'weight')
        parameters.weight = sum(molweight(fieldnames(r)));
        parameters.masses = molweight(fieldnames(r));
      end
      [b_coh, b_inc,sigma_abs,elements] = sqw_phonons_b_coh(fieldnames(r));
      % b     = sqrt(sigma*100/4/pi) [fm]
      % sigma = 4*pi*b^2/100         [barn]
      parameters.b_coh     = b_coh;    % [fm]
      parameters.b_inc     = b_inc;    % [fm]
      parameters.sigma_coh = sum(4*pi*b_coh.^2/100);
      parameters.sigma_inc = sum(4*pi*b_inc.^2/100);
      parameters.sigma_abs = sum(sigma_abs);
      parameters.atoms     = fieldnames(r);
      parameters.Concentration = cell2mat(struct2cell(r));
      parameters.At_number = zeros(size(parameters.atoms));
      for index=1:numel(parameters.atoms)
         z = find(strcmp(Z, parameters.atoms{index}),1);
         if ~isempty(z)
           parameters.At_number(index) = z;
         end
      end
    catch ME
      % disp([ mfilename ': invalid ChemicalFormula "' parameters.ChemicalFormula '" removed.' ]);
      parameters = rmfield(parameters, 'ChemicalFormula');
    end
  end
  
  if isfield(parameters, 'ChemicalFormula')
    s = setalias(s, 'ChemicalFormula', 'parameters.ChemicalFormula', ...
        'Material Chemical Formula');
  end
