function self=Bosify(self, type)
% iFunc_Sqw2D: Bosify: apply the 'Bose' factor (detailed balance) to a classical S(q,w) model.
%   The initial model should obey S*=S(q,w) = S(q,-w), i.e. be 'classical'.
%   The resulting model is 'quantum/experimental' and satisfies the detailed
%   balance. It contains the temperature effect (population).
%
% conventions:
% omega = Ei-Ef = energy lost by the neutron, given in [meV]
%    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
%    omega < 0, neutron gains energy, anti-Stokes
%
%    S(q,-w) = exp(-hw/kT) S(q,w)
%    S(q,w)  = exp( hw/kT) S(q,-w)
%    S(q,w)  = Q(w) S*(q,w) with S*=classical limit and Q(w) defined below.
%    for omega > 0, S(q,w) > S(q,-w)
%         
% The semi-classical correction, Q, aka 'quantum' correction factor, 
% can be selected from the optional   'type' argument:
%    Q = exp(hw/kT/2)                 'Schofield' or 'Boltzmann'
%    Q = hw/kT./(1-exp(-hw/kT))       'harmonic'  or 'Bader'
%    Q = 2./(1+exp(-hw/kT))           'standard'  or 'Frommhold' (default)
%
% The 'Boltzmann' correction leads to a divergence of the S(q,w) for energies above 
% kT (few 100 meV). The 'harmonic' correction provides a reasonable correction but
% does not fully avoid the divergence at large energies.
%
%  Bose factor: n(w) = 1./(exp(w*11.605/T) -1) ~ exp(-w*11.605/T)
%               w in [meV], T in [K]
%
% References:
%  B. Hehr, http://www.lib.ncsu.edu/resolver/1840.16/7422 PhD manuscript (2010).
%  S. A. Egorov, K. F. Everitt and J. L. Skinner. J. Phys. Chem., 103, 9494 (1999).
%  P. Schofield. Phys. Rev. Lett., 4, 239 (1960).
%  J. S. Bader and B. J. Berne. J. Chem. Phys., 100, 8359 (1994).
%  T. D. Hone and G. A. Voth. J. Chem. Phys., 121, 6412 (2004).
%  L. Frommhold. Collision-induced absorption in gases, 1 st ed., Cambridge
%    Monographs on Atomic, Molecular, and Chemical Physics, Vol. 2,
%    Cambridge Univ. Press: London (1993).
%
% syntax:
%   sqw_T = Bosify(sqw)
%   sqw_T = Bosify(sqw, type)
%
% input:
%   sqw:  Sqw model (classical, symmetric in energy, no T Bose factor)
%           e.g. 2D model with q as 1st axis (rows, Angs), w as 2nd axis (meV).
%   type: 'Schofield' or 'harmonic' or 'standard' (default)
%
% output:
%   sqw_T: quantum Sqw model (non classical, iFunc_Sqw2D).
%
% Example: s=sqw_difusion; sb=Bosify(s);
%
% See also: iFunc_Sqw2D, iData_Sqw2D
% (c) E.Farhi, ILL. License: EUPL.

  if nargin < 2, type = 'standard'; end
  
  % check if Temperature is already a parameter.
  T = findfield(self, 'Temperature','first');
  if isempty(T), T = findfield(self, 'T','exact first');
  % check that it is in the Parameters
  if ~isempty(T)  % re-use existing Temperature
    T_index = find(~cellfun(@isempty, strfind(self.Parameters, T)));
  else  % add new Temperature parameter
    self.Parameters{end+1} = 'Temperature [K]'; 
    T_index=numel(self.Parameters);
  end
  
  % check if the model is labelled as 'classical' or 'quantum'. get alias.
  classical = findfield(self, 'classical','first');
  if ~isempty(classical)
    classical = get(self, classical);
    if ~classical
      warning([ mfilename ': WARNING: Not "classical/symmetric": The model ' self.Tag ' ' self.Title ' does not seem to be classical (classical=0).' NL ...
  '   It may ALREADY contain the Bose factor in which case the detailed balance will be wrong.' ]);
    end
  end
  
  % handle deBosify as well
  if strfind(lower(type),'debosify')
    do_bosify=0;
    type = strtrim(strrep(lower(type), 'debosify','')); % remove debosify occurence
  else do_bosify=1; end
  
  % select the proper correction
  if strcmpi(type, 'Schofield') || strcmpi(type, 'exp') || strcmpi(type, 'Boltzmann') 
    % P. Schofield. Phys. Rev. Lett., 4, 239 (1960).
    Q  = 'exp(hw_kT/2)';          
  elseif strcmpi(type, 'Harmonic') || strcmpi(type, 'Bader')
    % J. S. Bader and B. J. Berne. J. Chem. Phys., 100, 8359 (1994).
    % p.t. D. Hone and G. A. Voth. J. Chem. Phys., 121, 6412 (2004).
    Q  = 'hw_kT./(1-exp(-hw_kT));  % w*(1+n(w))';
  elseif strcmpi(type, 'Standard') || strcmpi(type, 'Frommhold') || strcmpi(type, 'default')
    % L. Frommhold. Collision-induced absorption in gases, 1 st ed., Cambridge
    %   Monographs on Atomic, Molecular, and Chemical Physics, Vol. 2,
    %   Cambridge Univ. Press: London (1993).
    Q = '2./(1+exp(-hw_kT))';
  else
    error([ mfilename ': Bosify: Unknown semi-classical correction ' type ]);
  end
  
  % append the quantum correction to the Expression
  % apply sqrt(Bose) factor to get experimental-like
  % semi-classical corrections, aka quantum correction factor
  self.Expression{end+1} = [ 'kT = p(' num2str(T_index) ')/11.6045; % Kelvin to meV = 1000*K_B/e' ];               
  self.Expression{end+1} = 'hw_kT = y/kT; % hbar omega / kT';
  self.Expression{end+1} = [ 'Q=' Q ';' ];
  self.Expression{end+1} = 'Q(find(hw_kT==0)) = 1;';
  if ~do_bosify
    self.Expression{end+1} = 'Q = 1./Q;';
  end
  self.Expression{end+1} = 'signal = signal .* Q;  % apply detailed balance with the selected correction';
  
  if isvector(self.Guess) && isnumeric(self.Guess)
    self.Guess(T_index) = 300;
  end
  
  if nargout == 0 && ~isempty(inputname(1)) && isa(self,'iFunc')
    assignin('caller',inputname(1),self);
  end

end % Bosify
