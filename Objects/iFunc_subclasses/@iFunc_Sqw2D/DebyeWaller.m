function self=DebyeWaller(self)
% iFunc_Sqw2D: DebyeWaller: add the Debye-Waller factor to a S(q,w) model.
%   It accounts for thermal motions around the equilibrium.
%
%  The Debye-Waller factor for an harmonic, isotropic material is:
%    DW = exp(-u2*q^2/3)
%  where u2 is the mean squared displacement (in Angs^2). It gives the fraction 
%  of pure elastic scattering, whereas 1-DW is the pure inelastic scattering fraction.
%  The DW factor is related to the B factor or temperature factor:
%    B = 8 pi^2 <u2>
%  which is common in diffraction/Rietveld refinement.
%
% New model parameters:
%       u2     Mean squared displacement [Angs^2]
%
% References:
%  H. Schober, J. Neut. Res. 17 (2014) 109-357
%  G.L. Squires, Introduction to the Theory of Thermal Neutron Scattering, Dover Publications Inc.(1997)
%  <https://en.wikipedia.org/wiki/Debye%E2%80%93Waller_factor>
%
% syntax:
%   sDW = DebyeWaller(sqw)
%
% input:
%   sqw:  Sqw model
%           e.g. 2D model with q as 1st axis (rows, Angs), w as 2nd axis (meV).
%
% output:
%   sDW: Sqw model with DW factor (iFunc_Sqw2D).
%
% Example: s=sqw_difusion; DebyeWaller(s);
%
% See also: iFunc_Sqw2D, iData_Sqw2D
% (c) E.Farhi, ILL. License: EUPL.

  % search for 'msd','u2' parameters

  % add new <u2> parameter
  if isvector(self.Guess) && isnumeric(self.Guess)
    self.Guess(end+1) = 0.05;  % typical
  else
    if ~iscell(self.Guess), self.Guess = { self.Guess }; end
    self.Guess{end+1} = 0.05;
  end
  
  self.Parameters{end+1} = [ 'u2 Mean squared displacement [Angs^2] ' mfilename ]; 
  
  % check if the model has already a Debye-Waller factor ?
  hasDW = findfield(self, 'DebyeWaller','first');
  wrn = false;
  if ~isempty(hasDW)
    hasDW = get(self, hasDW);
    if hasDW(1)
      warning([ mfilename ': WARNING: The model ' self.Tag ' ' self.Name ' seems to already contain a Debye-Waller factor.' ]);
      wrn = true;
    end
  end
  
  t = self.Name;
  
  % append the quantum correction to the Expression
  % apply sqrt(Bose) factor to get experimental-like
  % semi-classical corrections, aka quantum correction factor
  self.Expression{end+1} = [ 'u2 = p(' num2str(numel(self.Parameters)) ');' ];               
  self.Expression{end+1} = 'DW = exp(-u2.*x.^2/3);';
  self.Expression{end+1} = 'signal = signal .* DW;';
  
  self.UserData.DebyeWaller     = true; 
  self.Name = [ mfilename '(' t ')' ];
  
  self = iFunc_Sqw2D(self); % check
  
  if nargout == 0 && ~isempty(inputname(1)) && isa(self,'iFunc') && ~wrn
    assignin('caller',inputname(1),self);
  end

end % DebyeWaller
