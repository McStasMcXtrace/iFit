function [t, fig]=thermochemistry(s, T, options)
% iFunc_Sqw4D/thermochemistry: compute thermodynamic quantities for 4D S(q,w) models.
%
% compute thermodynamics quantities for a 4D S(q,w) Model.
% When missing, the phonon (vibrational) density of states (vDOS) is estimated 
% from e.g. the evaluation of the 4D sqw phonons object.
%
% The input 4D model should be e.g. a 4D phonon one, with axes qh,qk,ql,energy
% in [rlu^3,meV]. The evaluation of the model should provide the density of states
% in model.UserData.DOS as an iData object.
%
% the function returns an array of iData objects:
%   DOS                                 [modes/energy unit/unit cell]
%   DOS_partials                    DOS partials (one per mode)
%   entropy                           S [eV/K/cell]
%   internal_energy                   U [eV/cell]
%   helmholtz_energy                  F [eV/cell]
%   heat_capacity at constant volume Cv [eV/K/cell]
%
% input:
%   s: a 4D sqw phonons model [K,H,L,Energy] (iFunc)
%   T: temperature range. When not given, T=1:500 [K]
%   options: can be 'plot' to display the results.
%         The 'newplot' option does the same but opens a new figure instead of
%         re-using an existing one.
%
% output:
%   t: structure with [DOS, DOS_partials, S U F Cv ]
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
if ~isa(s,'iFunc') || ndims(s) ~= 4
  disp([ mfilename ': Invalid model dimension. Should be iFunc/iData 2D or 4D. It is currently ' class(s) ' ' num2str(ndims(s)) 'D' ]);
  return
end

if isempty(T)
  T = 1:500;  % default
end

DOS = dos(s); % compute or get already available DOS from 4D

if nargout == 0 || ~isempty(strfind(options, 'plot'))
  options = [ options ' plot' ];
end

t = thermochemistry(DOS, T, options); % we use the iData_vDOS/thermochemistry

