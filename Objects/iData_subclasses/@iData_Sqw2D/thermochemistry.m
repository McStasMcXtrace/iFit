function [t, fig]=thermochemistry(s, T, options)
% iData_Sqw2D/thermochemistry: compute thermodynamic quantities for 2D S(q,w) data sets.
%
% compute thermodynamics quantities.
% When missing, the generalised density of states (gDOS) is estimated 
% from the 2D sqw object.
%
% The input 2D data set should be e.g. a 2D S(q,w), with axes |q|,energy
% in [Angs-1,meV]. 
%
% the function returns an array of iData objects:
%   DOS                                 [modes/energy unit/unit cell]
%   entropy                           S [eV/K/cell]
%   internal_energy                   U [eV/cell]
%   helmholtz_energy                  F [eV/cell]
%   heat_capacity at constant volume Cv [eV/K/cell]
%
% input:
%   s: a 2D sqw data set [|Q|,Energy] (iData_Sqw2D)
%   T: temperature range. When not given, T=1:500 [K]
%   options: can be 'plot' to display the results.
%         The 'newplot' option does the same but opens a new figure instead of
%         re-using an existing one.
%
% output:
%   t: structure with [DOS, S U F Cv ]
%   fig: figure handle
%
% Reference: https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html#background
%            D.A. McQuarrie. Statistical Mechanics. University Science Books, 2000.

t=[]; fig=[];
% test the input object
if nargin == 0, return; end
if nargin < 2, T=[]; end
if nargin < 3, options=[]; end

% must be 2D or 4D iFunc/iData.
if ~isa(s,'iData') || ndims(s) ~= 2
  disp([ mfilename ': Invalid model dimension. Should be iData 2D . It is currently ' class(s) ' ' num2str(ndims(s)) 'D' ]);
  return
end

if isempty(T)
  T = 1:500;  % default
end

DOS = dos(s);
omega_e = DOS{1};
dos_e   = DOS{0};

% renormalize the dos in eV
omega_e = omega_e/1000; % meV -> eV
dos_e   = dos_e*1000;

UD = s.UserData;

% factor to convert from eV/cell to J/mol
mole = 6.022140857e23;
e    = 1.60217733e-19;         % elementary charge,  C
factor = e*mole;               % [J/mol]
titl   = strtok(s.Name);
if isfield(s.UserData,'properties') && isfield(s.UserData.properties,'chemical_formula')
  chem = [ ' ' s.UserData.properties.chemical_formula ];
  titl = [ titl chem ];
else chem = ''; end
if isfield(s.UserData,'properties') && isfield(s.UserData.properties,'molar_mass')
  mass = sprintf('/%.1f [g]', s.UserData.properties.molar_mass);
  titl = [ titl mass ];
else mass = ''; end

% compute the entropy S [eV/K]
S = get_entropy(omega_e, dos_e, T)*factor;
entropy = iData(T, S);
entropy.Title=[ 'Entropy S=-dF/dT [J/K/mol] ' titl ]; 
xlabel(entropy,'Temperature [K]'); ylabel(entropy, [ 'S [J/K/mol]' chem mass ]);
entropy.Error=0; UD.entropy = entropy;

% compute the internal energy U [eV]
if isfield(s.UserData,'properties') && isfield(s.UserData.properties,'potential_energy')
  potential_energy = s.UserData.properties.potential_energy;
else potential_energy = 0; end

U = get_internal_energy(omega_e, dos_e, T, potential_energy)*factor;
internal_energy = iData(T, U);
internal_energy.Title=[ 'Internal energy U [J/mol] ' titl ]; 
xlabel(internal_energy,'Temperature [K]'); ylabel(internal_energy,[ 'U [J/moll]' chem mass ]);
internal_energy.Error=0; UD.internal_energy = internal_energy;

% compute the Helmotlz free energy F [eV]
F = U - T .* S;
helmholtz_energy = iData(T, F);
helmholtz_energy.Title=[ 'Helmholtz free energy F=U-TS=-kT lnZ [J/mol] ' titl ]; 
xlabel(helmholtz_energy,'Temperature [K]'); ylabel(helmholtz_energy,[ 'F [J/mol]' chem mass ]);
helmholtz_energy.Error=0; UD.helmholtz_energy = helmholtz_energy;

% compute the Cv
Cv = diff(internal_energy);
Cv.Title=[ 'Specific heat at constant volume Cv=dU/dT [J/K/mol] ' titl ]; 
xlabel(Cv,'Temperature [K]'); ylabel(Cv, [ 'Cv [J/K/mol]' chem mass ]);
Cv.Error=0; UD.heat_capacity = Cv;

s.UserData = UD;

if ~isempty(inputname(1))
  assignin('caller',inputname(1),s);
end

% return all
t.README = [ mfilename ' ' s.Name ];
t.DOS             =DOS;
t.Temperature     =T;
t.maxFreq         =s.UserData.maxFreq;
t.entropy         =s.UserData.entropy;            % S
t.internal_energy =s.UserData.internal_energy;    % U
t.helmholtz_energy=s.UserData.helmholtz_energy;   % F = U-TS
t.heat_capacity   =s.UserData.heat_capacity;      % Cv = dU/dT

if nargout == 0 || ~isempty(strfind(options, 'plot'))
  if ~isempty(strfind(options, 'newplot')) fig = figure; else fig = gcf; end
  subplot([ t.entropy t.internal_energy t.helmholtz_energy t.heat_capacity ], [3 2]);
  hold on
  if isfield(s.UserData,'DOS') && ~isempty(s.UserData.DOS)
      subplot(3,2,5);
      DOS = s.UserData.DOS;
      DOS{1} = DOS{1}; % change energy unit
      xlabel(DOS,[ 'Energy [meV]' ]);
      % plot any partials first
      
      % plot total DOS and rotate
      h=plot(DOS); set(h,'LineWidth',2);
    end
    set(fig, 'NextPlot','new');
end


% ------------------------------------------------------------------------------
function U = get_internal_energy(omega_e, dos_e, T, potential_energy)
  % get the internal energy U in eV/unit cell
  % 
  % omega_e: in eV
  
  if nargin < 3, T=[]; end
  if nargin < 4, potential_energy=0; end
  if isempty(T), T=300; end

  U = [];
  if numel(T) > 1
    for index=1:numel(T)
      U = [ U get_internal_energy(omega_e, dos_e, T(index), potential_energy) ];
    end
    return
  end
    
  U = potential_energy;
  % zero-point-energy
  zpe = trapz(omega_e, omega_e.*dos_e/2);

  k = 1.380658e-23;           % Boltzmann constant, J/K
  e = 1.60217733e-19;         % elementary charge,  C
  kB        = k / e;          % Boltzmann constant, eV/K
  B         = 1/(kB*T);       % beta=1/kB.T
    
  U         = U + zpe;
  E_vib     = omega_e ./ (exp(omega_e * B) - 1);
  E_vib(~isfinite(E_vib)) = 0;
  E_phonon  = trapz(omega_e, E_vib .* dos_e);
  U         = U + E_phonon;

function S = get_entropy(omega_e, dos_e, T)
  % get the entropy in eV/K/unit cell
  % 
  % omega_e: in eV
  
  if nargin < 3, T=[]; end
  if isempty(T), T=300; end
  
  S = [];
  if numel(T) > 1
    for index=1:numel(T)
      S = [ S get_entropy(omega_e, dos_e, T(index)) ];
    end
    return
  end
  
  k = 1.380658e-23;           % Boltzmann constant, J/K
  e = 1.60217733e-19;         % elementary charge,  C
  kB        = k / e;          % Boltzmann constant, eV/K
  B         = 1/(kB*T);       % beta=1/kB.T
  
  S_vib = omega_e ./ (T * (exp(omega_e * B) - 1)) - kB * log(1 - exp(-omega_e * B));
  S_vib(~isfinite(S_vib)) = 0;
                 
  S = trapz(omega_e, S_vib .* dos_e);
  

