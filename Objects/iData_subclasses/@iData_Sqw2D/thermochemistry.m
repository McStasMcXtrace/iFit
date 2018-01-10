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

DOS = dos(s); % compute or get already available DOS from 2D

if nargout == 0 || ~isempty(strfind(options, 'plot'))
  options = [ options ' plot' ];
end

t = thermochemistry(DOS, T, options); % we use the iData_vDOS/thermochemistry
