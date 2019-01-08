function signal=sqw_visco_elastic(varargin)
% model = sqw_visco_elastic(p, q, w, {signal}) : Visco Elastic model
%
%  iFunc/sqw_visco_elastic Visco Elastic model in the linear response theory of liquids.
%    This model can be used to fit e.g. a spectrum with a central non dispersive 
%    quasi elastic line  (lorentzian shape) and two dispersive (e.g. acoustic)
%    Stokes and anti-Stokes lines. This model is suited for e.g dense fluids at 
%    small Q momentum values.
%  The expression is derived from the Egelstaff longitudinal modes theory in 
%   liquids (Eq 15.27).
%
%  S(q,w)/S(q) = (1-1/CpCv) * Dq2 / (Dq2^2 + w^2) 
%              + 1/2/CpCv*gq2* (1./((w+cq)^2+gq2^2) + 1/((w-cq)^2+gq2^2) )
%
%  The diffusion constant D is usually around D=1-10 E-9 [m^2/s] in liquids. Its 
%  value in [meV/Angs^2] is D*4.1356e+08.The sound velocity is usually around 
%  1000 [m/s]. Its value in [meV/Angs] is c/142.1622.
%
%  In this implementation, S(q) is assumed to be 1, however, it is possible to 
%  set S(q) over a full q-range when creating the model:
%    s = sqw_visco_elastic(sq)
%  where 'sq' is an iData [q,s(q)] in [Angs^-1], which is stored in UserData.Sq
%  The S(q) can be also changed or set anytime as:
%    s.UserData.Sq = iData(q,Sq);
%  It is also possible to multiply this model with a S(q) model (see Examples below).
%
%  The S(q,w) is computed in its symmetrized expression (so-called classical). 
%  To get the 'true' quantum S(q,w) definition, use e.g.
%    sqw = sqw_visco_elastic .* bose';
%  where the Temperature is then given in [x units]. If 'x' is an energy in [meV]
%  then the Temperature parameter is T[K]/11.6045
%
%  To add a 'background' use e.g.
%    sqw = sqw_visco_elastic + constant('Background');
%
%  To convolute with a pseudo-Voigt or Gaussian instrumental line width, use e.g.
%    sqw = convn(sqw_visco_elastic, gauss2d)
% or
%    sqw = convn(sqw_visco_elastic, pseudovoigt2d)
%
% Reference: 
%  P.A.Egelstaff, An introduction to the liquid state, 2nd ed., Oxford (2002)
%
% Example:simple liquid with d=1 Angs inter-atomic distance and static S(q)
%  sq = iData(sf_hard_spheres, [1 0.4], 0:0.01:20); % a simple liquid with d=1 Angs
%  s  = sqw_visco_elastic(sq);
%  plot(log10(iData(s, [], 0:.01:4, -50:50)))
%
% Example: simple liquid with parameterized S(q)
%  s1 = sqw_visco_elastic; s1.ParameterValues=s1.Guess;
%  s2 = sf_hard_spheres;   s2.ParameterValues=[4 .4];
%  s  = s1.*s2; 
%  plot(log10(iData(s, [], 0:.01:4, -50:50)))
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot, gauss, sqw_phonons, sqw_cubic_monoatomic, sqw_vaks
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>, sqw_gen_hydrodynamics
% (c) E.Farhi, ILL. License: EUPL.

signal.Name           = [ 'sqw_visco_elastic Visco Elastic - Rayleight-Brillouin triplet [' mfilename ']' ];
signal.Description    = 'Visco Elastic model with a central line, and two dispersive acoustic lines, suited for e.g. dense fluids';

signal.Parameters     = {  ...
  'Amplitude' ...
  'c      sound velocity in [m/s]' ...
  'Dt     Diffusion coefficient for thermal central line [m^2/s]' ...
  'Dl     Diffusion coefficient for longitudinal mode [m^2/s]' ...
  'gamma  Laplace coefficient Cp/Cv [1]' ...
};
  
signal.Dimension      = 2;         % dimensionality of input space (axes) and result
signal.UserData.classical     = true;
signal.UserData.DebyeWaller   = false;

% we must provide a p(xx) expression to tell iFunc to validate xx parameters
signal.Expression = { ...
  '% the energy axis' ...
  'q = x; w=y; c=p(2); D=p(3); g=p(4); CpCv=p(5);' ...
  'if isvector(q) && isvector(w) && numel(q) ~= numel(w), [q,w] = meshgrid(q,w); end' ...
  'dq2 = D*q.^2*1E20/241.8E9; % in meV, with q in Angs-1 and D in m2/s' ...
  'cq  = c*q/142.1622;        % m/s -> meV/Angs' ...
  'gq2 = g*q.^2*1E20/241.8E9; % in meV, with q in Angs-1 and D in m2/s' ...
  'signal = (1-1/CpCv) * dq2 ./ (dq2.^2 + w.^2);' ...
  'signal = signal + 1/2/CpCv.*gq2.* (1./((w+cq).^2+gq2.^2) + 1./((w-cq).^2+gq2.^2) );' ...  
  'signal = signal*p(1)/pi; ' ...
  'if isfield(this.UserData,''Sq'')' ...
    'sq = this.UserData.Sq;' ...
    'if isa(sq, ''iData'') && ndims(sq) == 1' ...
      'try' ...
      'sq = double(interp(sq, unique(x(:)))); sq = transpose(sq(:));' ...
      'sq = repmat(sq, numel(unique(y(:))),1);' ...
      'signal = signal .* sq;' ...
      'end' ...
    'end' ...
  'end' ...
};

signal.Guess          = [ 1 500 1e-9 1e-9 1.5 ]; % default parameters

signal=iFunc(signal);
signal=iFunc_Sqw2D(signal);

% treat input argument at model creation
for index=1:nargin
  varg = varargin{index};
  if isa(varg, 'iData') && ndims(varg) == 1, signal.UserData.Sq = varg;
  elseif isnumeric(varg) && isvector(varg)
    signal.ParameterValues = varg;
  end
end

