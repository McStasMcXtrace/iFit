classdef iData_vDOS < iData

  % iData_vDOS: create a vibrational density of states data set (iData flavour)
  %
  % The iData_vDOS class is a 1D data set holding a vibrational density of states
  % often labelled as g(w)
  %
  % Example: g=iData_vDOS('SQW_coh_lGe.nc')
  %
  % Useful methods for this iData flavour:
  %
  % methods(iData_Sqw2D)
  %   all iData methods can be used.
  % iData_vDOS(g)
  %   convert input [e.g. a 2D or 4D Sqw, or a 1D vDOS] into an iData_vDOS to give access to
  %   the methods below.
  % d = dos(g)
  %   Compute/get the vibrational density of states (vDOS).
  % t = thermochemistry(g)
  %   Compute and display thermochemistry quantities from the vDOS.
  % [Sqw,Iqt,DW] = multi_phonons_incoherent(g, q, T, sigma, m, n)
  %   Compute the so-called multi-phonon dynamic structure factor expansion in 
  %   the incoherent gaussian approximation. The neutron scattering incoherent 
  %   gaussian dynamic structure factor is computed, as well as the corresponding
  %   intermediate scattering function, and the Debye-Waller factor.
  % [Gw, Tsym] = multi_phonons_vDOS(g, Ki, T, sigma, m, phi, n)
  %   Compute the so-called multi-phonon density of states expansion in the 
  %   neutron scattering incoherent gaussian approximation.
  %
  % input:
  %   can be an iData (1D vDOS, 2D Sqw, 4D Sqw) or filename or 2D/4D Sqw Model
  %   to generate a 1D vDOS object.
  %
  % output: an iData_vDOS object
  %
  % See also: iData

  properties
  end
  
  methods
  end
  
end
