% Calculate sqw file header and data for a single spe file
%
%   >> [header,sqw_data]=calc_sqw (efix, emode, alatt, angdeg, u, v,...
%                                       psi, omega, dpsi, gl, gs, data, det)
% Input:
%   efix        Fixed energy (meV)
%   emode       Direct geometry=1, indirect geometry=2
%   alatt       Lattice parameters (Ang^-1)
%   angdeg      Lattice angles (deg)
%   u           First vector (1x3) defining scattering plane (r.l.u.)
%   v           Second vector (1x3) defining scattering plane (r.l.u.)
%   psi         Angle of u w.r.t. ki (rad)
%   omega       Angle of axis of small goniometer arc w.r.t. notional u
%   dpsi        Correction to psi (rad)
%   gl          Large goniometer arc angle (rad)
%   gs          Small goniometer arc angle (rad)
%   data        Data structure of spe file (see get_spe)
%            or The same, but with in addition a field qspec, a 4xn array of qx,qy,qz,eps
%   det         Data structure of par file (see get_par)
%            or If data has field qspec, this is ignored
%
% Ouput:
%   header      Header information in data structure suitable for write_sqw_header
%   sqw_data    Data structure suitable for write_sqw_data
%