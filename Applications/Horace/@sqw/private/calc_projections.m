% Label pixels in an spe file with coords in the 4D space defined by crystal Cartesian coordinates
% and energy transfer. 
% Allows for correction scattering plane (omega, dpsi, gl, gs) - see Tobyfit for conventions
%
%   >> [u_to_rlu, ucoords] = ...
%    calc_projections (efix, emode, alatt, angdeg, u, v, psi, omega, dpsi, gl, gs, data, det)
%
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
% Output:
%   u_to_rlu    Matrix (3x3) of projection axes in reciprocal lattice units
%              i.e. u(:,1) first vector - u(1:3,1) r.l.u. etc.
%              This matrix can be used to convert components of a vector in the
%              projection axes to r.l.u.: v_rlu = u * v_proj
%   ucoords     [4 x npix] array of coordinates of pixels in crystal Cartesian
%              coordinates and energy transfer
%