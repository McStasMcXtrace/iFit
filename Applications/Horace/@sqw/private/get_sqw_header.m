% Read the header block for the results of performing calculate projections on spe file(s).
%
% Syntax:
%   >> [data, mess] = get_sqw_header(fid, data_in)
%
% Input:
% ------
%   fid         File pointer to (already open) binary file
%
% Output:
% -------
%   data        Structure containing fields read from file (details below)
%   mess        Error message; blank if no errors, non-blank otherwise
%
% Fields read from file are:
%   data.filename   Name of sqw file excluding path
%   data.filepath   Path to sqw file including terminating file separator
%   data.efix       Fixed energy (ei or ef depending on emode)
%   data.emode      Emode=1 direct geometry, =2 indirect geometry
%   data.alatt      Lattice parameters (Angstroms)
%   data.angdeg     Lattice angles (deg)
%   data.cu         First vector defining scattering plane (r.l.u.)
%   data.cv         Second vector defining scattering plane (r.l.u.)
%   data.psi        Orientation angle (deg)
%   data.omega      --|
%   data.dpsi         |  Crystal misorientation description (deg)
%   data.gl           |  (See notes elsewhere e.g. Tobyfit manual
%   data.gs         --|
%   data.en         Energy bin boundaries (meV) [column vector]
%   data.uoffset    Offset of origin of projection axes in r.l.u. and energy ie. [h; k; l; en] [column vector]
%   data.u_to_rlu   Matrix (4x4) of projection axes in hkle representation
%                      u(:,1) first vector - u(1:3,1) r.l.u., u(4,1) energy etc.
%   data.ulen       Length of projection axes vectors in Ang^-1 or meV [row vector]
%   data.ulabel     Labels of the projection axes [1x4 cell array of character strings]
%