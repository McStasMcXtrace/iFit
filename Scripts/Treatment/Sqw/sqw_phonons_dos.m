function [s,f]=sqw_phonons_dos(s, T)
% sqw_phonons_dos: compute the vDOS and thermodynamic quantities
%
% compute the vDOS from the evaluation of the 4D sqw_phonons object
% then compute thermodynamics quantities.
% The input 4D model should be e.g. a 4D phonon one, with axes qh,qk,ql,energy
% in [rlu^3,meV]. The evaluation of he model should provide the bare excitation
% frequencies into the UserData.FREQ, in [meV].
%
% input:
%   s: a 4D sqw_phonons model (iFunc)
%
% output:
%   s: enriched model
%   f: iData holding the 

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
  qh=linspace(0.01,.5,10);qk=qh; ql=qh; w=linspace(0.01,50,11);
  f=iData(s,[],qh,qk,ql',w);
  s.UserData.maxFreq = max(s.UserData.FREQ(:));
  disp([ mfilename ': maximum phonon energy ' num2str(s.UserData.maxFreq) ' [meV] in ' s.Name ]);
end 

% re-evaluate the 4D model onto a mesh filling the Brillouin zone [-0.5:0.5]
s.UserData.maxFreq = max(s.UserData.maxFreq(:));
qh=linspace(-0.5,.5,30);qk=qh; ql=qh; w=linspace(0.01,s.UserData.maxFreq*1.2,51);
f=iData(s,[],qh,qk,ql',w);

% then get back the bare frequencies
if ~isfield(s.UserData,'FREQ')
  disp([ mfilename ': The object should have the bare frequencies available as UserData.FREQ [meV].' ]);
  return
end
FREQ = s.UserData.FREQ; % in meV

% compute the DOS
index= find(FREQ > 0);
[dos_e,omega_e]=hist(FREQ(index),50);
Escaling = findfield(s, 'Energy_scaling', 'first');
if ~isempty(Escaling), Escaling = get(s, Escaling); end
if isnumeric(Escaling) && Escaling > 0, omega_e=omega_e*Escaling(1); end
DOS=iData(omega_e,dos_e./sum(dos_e));
DOS.Title = [ 'Total DOS ' strtok(s.Name) ]; 
xlabel(DOS,'Energy [meV]'); 
ylabel(DOS,'Total DOS');
UD = s.UserData;
DOS.Error=0; UD.DOS=DOS;

% clean up the energy axis
omega_e = omega_e/1000; % eV

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

U = get_internal_energy(omega_e, dos_e, T, potential_energy);
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

% update the object
if nargout == 0 && ~isempty(inputname(1))
  assignin('caller',inputname(1),s);
end

function D = get_total_dos(s, npts, delta)
  % compute the total Phonon Density of States
  if nargin < 2, npts = []; end
  if nargin < 3, delta= []; end
  if isempty(npts), npts=1000; end
  if isempty(delta), delta = 0.5; end

  omega_kl = s.UserData.FREQ; % in meV
  omega_e  = linspace(0,max(FREQ(:))*1.2, npts);
  dos_e    = zeros(size(omega_e));
  delta    = 0.05; % meV smoothing

  omega_el = repmat(omega_e, [3 1])';

  for index=1:size(omega_kl,1)  % along q values
    diff_el = ( omega_el - repmat(omega_kl(index,:),[npts 1]) ).^2;
    dos_el  = 1./ (diff_el + (0.5*delta^2)); 
    dos_e   = dos_e + sum(dos_el,2)';
  end

% ------------------------------------------------------------------------------
function U = get_internal_energy(omega_e, dos_e, T, potential_energy)
  % get the internal energy U in eV/unit cell
  
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
  zpe = trapz(omega_e.*dos_e/2, omega_e);

  T2E       = (1/11.6045);           % Kelvin to meV = 1000*K_B/e
  kT        = T*T2E;
  hw_kT     = omega_e./kT;           % hbar omega / kT
    
  U         = U + zpe;
  E_vib     = omega_e ./ (exp(hw_kT) - 1);
  E_phonon  = trapz(E_vib .* dos_e, omega_e);
  U         = U + E_phonon;

function S = get_entropy(omega_e, dos_e, T)
  % get the entropy in eV/K/unit cell
  
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
  kT        = T*T2E;
  hw_kT     = omega_e./kT;    % hbar omega / kT
  k = 1.380658e-23;           % Boltzmann constant, J/K
  e = 1.60217733e-19;         % elementary charge,  C
  kB        = k / e;          % Boltzmann constant, eV/K
  
  S_vib = (omega_e ./ (T * (exp(hw_kT) - 1.)) - kB * log(1. - exp(-hw_kT)));
                 
  S = trapz(S_vib .* dos_e, omega_e);

