% Calculate the components of Q in reference frame fixed w.r.t. spectrometer
%
%   >> qspec = calc_qspec (efix, k_to_e, emode, data, det)
%
% Input:
% ------
%   efix    Fixed energy (meV)
%   k_to_e  Constant in the relation energy (meV) = k_to_e *(wavevector^2)
%   emode   Direct geometry=1, indirect geometry=2, elastic=0
%           If elastic, then interprets energy bins as logarithm of wavelength (Ang)
%   data    Data structure of spe file (see get_spe)
%   det     Data structure of par file (see get_par)
%
% Output:
% -------
%   qspec(4,ne*ndet)    Momentum and energy transfer in spectrometer coordinates
%
%  Note: We sometimes use this routine with the energy bin boundaries replaced with 
%        bin centres i.e. have fudged the array data.en
%