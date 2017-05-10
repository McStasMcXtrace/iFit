function [s,f]=sqw_thermochemistry(s, T)
% sqw_thermochemistry: compute the vDOS and thermodynamic quantities.
%
% Compute the phonon (vibrational) density of states (vDOS) from e.g. the evaluation
% of the 4D sqw_phonons object. Then compute thermodynamics quantities.
%
% The input 4D model should be e.g. a 4D phonon one, with axes qh,qk,ql,energy
% in [rlu^3,meV]. The evaluation of the model should provide the density of states
% in model.UserData.DOS as an iData object.
%
% the function returns an array of iData objects:
%   DOS 
%   pDOS (one per mode)
%   entropy                           S [eV/K/cell]
%   internal_energy                   U [eV/cell]
%   helmholtz_energy                  F [eV/cell]
%   heat_capacity at constant volume Cv [eV/K/cell]
%
% input:
%   s: a 4D sqw_phonons model [K,H,L,Energy] (iFunc)
%   T: temperature range. When not given, T=1:500 [K]
%
% output:
%   s: [DOS, partial_DOS, S U F Cv ]
%   f: iData holding the 4D evaluation of the model.

f = []; 
% test the input object
if nargin == 0, return; end

if (isa(s,'iFunc') || isa(s,'iData')) && ndims(s) == 2
  disp([ 'Using Sqw_gDOS for ' class(s) ' ' num2str(ndims(s)) 'D objects.' ]);
  [gDOS, g] = Sqw_gDOS(s);
  UD = s.UserData;
  UD.DOS = gDOS; f=g;
  s.UserData = UD;
  return
end
  
  
if ndims(s) ~= 4 || ~isa(s,'iFunc')
  disp([ mfilename ': Invalid model dimension. Should be iFunc 4D. It is currently ' class(s) ' ' num2str(ndims(s)) 'D' ]);
  return
end

% first get a quick estimate of the max frequency
if ~isfield(s.UserData,'maxFreq')
  qh=linspace(-.5,.5,10);qk=qh; ql=qh; w=linspace(0.01,50,11);
  f=iData(s,[],qh,qk,ql',w);
  s.UserData.maxFreq = max(s.UserData.FREQ(:));
  disp([ mfilename ': maximum phonon energy ' num2str(s.UserData.maxFreq) ' [meV] in ' s.Name ]);
end 

% evaluate the 4D model onto a mesh filling the Brillouin zone [-0.5:0.5 ]
s.UserData.DOS     = [];  % make sure we re-evaluate again on a finer grid
s.UserData.maxFreq = max(s.UserData.maxFreq(:));
qh=linspace(-0.5,.5,30);qk=qh; ql=qh; w=linspace(0.01,s.UserData.maxFreq*1.2,51);
f=iData(s,[],qh,qk,ql',w);

% then get back the bare frequencies
if ~isfield(s.UserData,'DOS')
  disp([ mfilename ': The object has no density of states available. Can not compute thermo-chemistry.' ]);
  return
end

% get the DOS
DOS     = s.UserData.DOS;
omega_e = DOS{1};
dos_e   = DOS{0};

% renormalize the dos in eV
omega_e = omega_e/1000; % meV -> eV
dos_e   = dos_e*1000;

% get the temperature
if nargin < 2
  T = 1:500;  % default
  try
    T = Sqw_getT(s);  % a value is defined in the Parameters ?
  end
end

% compute the entropy S [eV/K]
S = get_entropy(omega_e, dos_e, T);
entropy = iData(T, S);
entropy.Title=[ 'Entropy S=-dF/dT [eV/K/cell] ' strtok(s.Name) ]; 
xlabel(entropy,'Temperature [K]'); ylabel(entropy, 'S [eV/K/cell]');
entropy.Error=0; UD.entropy = entropy;

% compute the internal energy U [eV]
if isfield(s.UserData,'properties') && isfield(s.UserData.properties,'potential_energy')
  potential_energy = s.UserData.properties.potential_energy;
else potential_energy = 0; end

U = get_internal_energy(omega_e/1000, dos_e, T, potential_energy);
internal_energy = iData(T, U);
internal_energy.Title=[ 'Internal energy U [eV/cell] ' strtok(s.Name) ]; 
xlabel(internal_energy,'Temperature [K]'); ylabel(internal_energy,'U [eV/cell]');
internal_energy.Error=0; UD.internal_energy = internal_energy;

% compute the Helmotlz free energy F [eV]
F = U - T .* S;
helmholtz_energy = iData(T, F);
helmholtz_energy.Title=[ 'Helmholtz free energy F=U-TS=-kT lnZ [eV/cell] ' strtok(s.Name) ]; 
xlabel(helmholtz_energy,'Temperature [K]'); ylabel(helmholtz_energy,'F [eV/cell]');
helmholtz_energy.Error=0; UD.helmholtz_energy = helmholtz_energy;

% compute the Cv
Cv = diff(internal_energy);
Cv.Title=[ 'Specific heat at constant volume Cv=dU/dT [eV/K/cell] ' strtok(s.Name) ]; 
xlabel(Cv,'Temperature [K]'); ylabel(Cv, 'Cv [eV/K/cell]');
Cv.Error=0; UD.heat_capacity = Cv;

s.UserData = UD;

if ~isempty(inputname(1))
  assignin('caller',inputname(1),s);
end

% return all
s = [s.UserData.DOS s.UserData.pDOS ...
    s.UserData.entropy s.UserData.internal_energy s.UserData.helmholtz_energy ...
    s.UserData.heat_capacity ];

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

  T2E       = (1/11.6045);           % Kelvin to meV = 1000*K_B/e
  kT        = T*T2E/1000;
  hw_kT     = omega_e./kT;           % hbar omega / kT
    
  U         = U + zpe;
  E_vib     = omega_e ./ (exp(hw_kT) - 1);
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
  
  T2E       = (1/11.6045);    % Kelvin to meV = 1000*K_B/e
  kT        = T*T2E/1000;
  hw_kT     = omega_e./kT;    % hbar omega / kT
  k = 1.380658e-23;           % Boltzmann constant, J/K
  e = 1.60217733e-19;         % elementary charge,  C
  kB        = k / e;          % Boltzmann constant, eV/K
  
  S_vib = (omega_e ./ (T * (exp(hw_kT) - 1.)) - kB * log(1. - exp(-hw_kT)));
                 
  S = trapz(omega_e, S_vib .* dos_e);

