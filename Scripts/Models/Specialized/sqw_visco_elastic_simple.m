function signal=sqw_visco_elastic_simple(varargin)
% model = sqw_visco_elastic_simple(p, q, w, {signal}) : Visco Elastic model for Lennard-Jones model Sqw2D
%
%  iFunc/sqw_visco_elastic_simple: S(q,w) Visco Elastic model in the linear response 
%    theory of liquids for a Lennard-Jones/Percus-Yevick simple model.
%    This model can be used to fit e.g. a spectrum with a central non dispersive 
%    quasi elastic line  (lorentzian shape) and two linear dispersive (e.g. acoustic)
%    Stokes and anti-Stokes lines. This model is suited for e.g dense fluids 
%    for coherent processes.
%  The expression is derived from the Egelstaff longitudinal modes theory in 
%    liquids (Eq 15.27) and a logitudinal dispersion from 
%
%  Model and parameters:
%  ---------------------
%
%  S(q,w)/S(q) = (1-1/CpCv) Dt q2 / ((Dt q2)^2 + w^2) 
%               + 1/2/CpCv  Dl q2 (1/((w+cq)^2+(Dl q2)^2) + 1/((w-cq)^2+(Dl q2)^2) )
%
%  where:
%    CpCv=gamma is the Laplace coefficient usually in 1.1-1.67
%    cq         is the acoustic linear longitudinal dispersion so that
%                 cq^2 = c0^2*q^2+w0^2*(1-3*sin(qR)/(qR) - 6*cos(qR)/(qR)^2+6*sin(qR)/(qR)^3 )
%               To suppress the acoustic mode fix c=0 or Dl=0.
%    Dt         is the entropy fluctuation (thermal) diffusion coefficient
%                 the central line half width is then Dt*q^2. To suppress it fix Dt=0
%    Dl         is the longitudinal diffusion coefficient
%                 the acoustic mode half width is then Dl*q^2. To suppress it fix Dl=0
%
%  Usually Dt = lambda/Cp.rho with lambda the thermal conductivity [W/m/K] and 
%  rho the material density [g/cm3]. Also, the sound velocity c is 1/sqrt(rho Chi)
%  where Chi is the isothermal compressibility [1/Pa].
%
%  The diffusion constant D is usually around D=1-10 E-9 [m^2/s] in liquids. Its 
%  value in [meV.Angs^2] is D*4.1356e+08.The sound velocity is usually around 
%  1000 [m/s]. Its value in [meV.Angs] is c/142.1622.
%
%  Structure factor:
%  -----------------
%
%  In this implementation, the structure factor is modelled from a Percus-Yevick
%  model for hard spheres (Lennard-Jones). The parameter 'eta' ranges from 0 
%  [ideal liquid S(q)=1] to 1 [pure crystalline solid with sharp Bragg peaks].
%
%  Additional remarks:
%  -------------------
%
%  The S(q,w) is computed in its symmetrized expression (so-called classical)
%    and does NOT satisfy the detailed balance.
%
%  To get the 'true' quantum S(q,w) definition, use e.g.
%
%    sqw = sqw_visco_elastic_simple .* bose;
%  where the Temperature is then given in [x units]. If 'x' is an energy in [meV]
%  then the Temperature parameter is T[K]/11.6045
%
%  To add a 'background' use e.g.
%    sqw = sqw_visco_elastic_simple + constant('Background');
%
%  To convolute with a pseudo-Voigt or Gaussian instrumental line width, use e.g.
%    sqw = convn(sqw_visco_elastic_simple, gauss2d)
% or
%    sqw = convn(sqw_visco_elastic_simple, pseudovoigt2d)
%
% Usage:
% ------
% s = sqw_visco_elastic_simple; 
% value = s(p,q,w); or value = iData(s,p,q,w)
%
% input:  p: sqw_visco_elastic_simple model parameters (double)
%             p(1)=Amplitude
%             p(2)=c      Sound velocity in [m/s]
%             p(3)=Dt     Diffusion coefficient for thermal central line [m^2/s]
%             p(4)=Dl     Diffusion coefficient for longitudinal mode [m^2/s]
%             p(5)=gamma  Laplace coefficient Cp/Cv [1]
%             p(6)=w0     Maximum acoustic phonon energy [meV]
%             p(7)=R      Inter-atomic/molecule distance [Angs]
%             p(8)=eta    Reduced density [0-1]
%         q:  axis along wavevector/momentum (row,double) [Angs-1]
%         w:  axis along energy (column,double) [meV]
% output: signal: model value [iFunc_Sqw2D]
%
% Reference: 
%  P.A.Egelstaff, An introduction to the liquid state, 2nd ed., Oxford (2002)
%  Copley and Lovesey, Rep. Prog. Phys. 38 (1975) 461
%
% Example:simple liquid
%  s  = sqw_visco_elastic_simple;
%  plot(log10(iData(s, [], 0:.01:4, -50:50)))
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot, gauss, sqw_phonons, sqw_cubic_monoatomic, sqw_vaks
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>, sqw_gen_hydrodynamics, sqw_visco_elastic
% (c) E.Farhi, ILL. License: EUPL.


% NOTES ************************************************************************
% Copley Lovesey Rep Prog Phys Eq 2.39 p476: w = f(q) acoustic dispersion/characteristic
%
% wl^2 = 3*q.^2*kT/m+w0^2*( 1-3*sin(q*R)./(q*R) - 6*cos(q*R)./(q*R).^2+6*sin(q*R)./(q*R).^3 )
%      ~ 3*q.^2*kT/m+w0^2.*(q*R).^2/4
%
% so c ~ sqrt(3*kT/m+w0^2*R^2/4)
%    c^2    = 3*kT/m+w0^2*R^2/4
%
% R = 3.4 Angs in l-Rb
%
% cq = sqrt(wl^2)
% Percus Yevick Copley Lovesey Eq 2.30 p 472. eta is unitless.
% Egelstaff: eta = pi rho sigma^3/6 Eq 5.50    rho in 1/Angs^3  rho.sigma^3 ~ 0.8
%   sigma=LJ finite distance at which the inter-particle potential is zero 
%   -> interatomic equilibrium distance 
% S(0)  = 1/alpha
% alpha = (1+2*eta)^2/(1-eta)^4
% Isotropic_Sqw line 1000 -> S(0): ChiT= S(0)/(kT*rho*1e30) [Pa-1] (Egelstaff  p201 Eq 10.21)

% treat input argument at model creation
for index=1:nargin
  varg = varargin{index};
  if isnumeric(varg) && isvector(varg)
    signal.ParameterValues = varg;
  end
end

signal.Name           = [ 'sqw_visco_elastic_simple Visco Elastic/Lennard-Jones - Rayleight-Brillouin triplet [' mfilename ']' ];
signal.Description    = 'Visco Elastic/Lennard-Jones model with a central line, and two dispersive acoustic lines, suited for e.g. dense fluids';

signal.Parameters     = {  ...
  'Amplitude' ...
  'c      Sound velocity in [m/s]' ...
  'Dt     Diffusion coefficient for thermal central line [m^2/s]' ...
  'Dl     Diffusion coefficient for longitudinal mode [m^2/s]' ...
  'gamma  Laplace coefficient Cp/Cv [1]' ...
  'w0     Maximum acoustic phonon energy [meV]' ...
  'R      Inter-atomic/molecule distance [Angs]' ...
  'eta    Reduced density [0-1]' ...
};
  
signal.Dimension      = 2;         % dimensionality of input space (axes) and result
signal.UserData.classical     = true;
signal.UserData.DebyeWaller   = false;

% we must provide a p(xx) expression to tell iFunc to validate xx parameters
signal.Expression = { ...
  '% the energy axis' ...
  'q = x; w=y; c=p(2); D=p(3); g=p(4); CpCv=p(5); w0=p(6); R=p(7); fp=p(8);' ...
  'if isvector(q) && isvector(w) && numel(q) ~= numel(w), [q,w] = meshgrid(q,w); end' ...
  'dq2 = D*q.^2*1E20/241.8E9; % in meV, with q in Angs-1 and D in m2/s' ...
  '% compute c0 constant c0^2=3*kT/m and c^2 = c0^2+w0^2*R^2/4' ...
  'c02 = (c/142.1622)^2 - w0^2*R^2/4; % m/s -> [meV.Angs]^2' ...
  'if c02 < 0, c02 = 0; end' ...
  'wl2 = c02*q.^2+w0^2*( 1-3*sin(q*R)./(q*R) - 6*cos(q*R)./(q*R).^2+6*sin(q*R)./(q*R).^3 );' ...
  'cq  = sqrt(wl2);' ...
  'gq2 = g*q.^2*1E20/241.8E9; % in meV, with q in Angs-1 and D in m2/s' ...
  'signal = 0;' ...
  'if D>0' ...
    'signal = (1-1/CpCv) * dq2 ./ (dq2.^2 + w.^2);' ...
  'end' ...
  'if g>0 && c>0' ...
    'signal = signal + 1/2/CpCv.*gq2.* (1./((w+cq).^2+gq2.^2) + 1./((w-cq).^2+gq2.^2) );' ...  
  'end' ...
  'signal = signal*p(1)/pi; ' ...
  '% use LJ S(q). Code from sf_hard_spheres' ...
  'fp=max(0, min(fp, 1));' ...
  'if (fp <= 0.0 || fp >= 1) sq=ones(size(q)); else ' ...
    'A=abs(R*q); '...
    'alpha = power(1.0+2.0*fp,2.0)/power(1.0-fp,4.0); ' ...
    'beta  = -6.0*fp*power(1.0+fp/2.0,2.0)/power(1.0-fp,4.0); gamma = fp*alpha/2.0; ' ...
    'sq = alpha*(sin(A)-A.*cos(A))./power(A,2.0); ' ...
    'sq = sq + beta*(2.0*A.*sin(A)+(2.0-power(A,2.0)) .* cos(A)-2.0)./power(A,3.0); ' ...
    'sq = sq + gamma * (-power(A,4.0) .* cos(A) + 4.0*((3.0*power(A,2.0)-6.0).*cos(A)+(power(A,3.0)-6.0*A).*sin(A)+6.0))./power(A,5.0); ' ...
    'sq = 1.0./(1.0+24.0*fp*sq./A); sq(isnan(sq))=0; end' ...
  'signal = signal .* sq;' ...
};

signal.Guess          = [ 1 2500 1e-9 1e-9 1.5 3.4 10 0.4 ]; % default parameters

signal=iFunc(signal);
signal=iFunc_Sqw2D(signal);



