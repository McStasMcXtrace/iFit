% Read one or more spe files and a detector parameter file, and create an output sqw file.
%
% Normal use:
%   >> accumulate_sqw (dummy, spe_file, par_file, sqw_file, efix, emode, alatt, angdeg,...
%                                               u, v, psi, psi_planned, omega, dpsi, gl, gs, grid_size_in, urange_in)
% If want output diagnostics:
%   >> [tmp_file,grid_size,urange] = gen_sqw (dummy, spe_file, par_file, sqw_file, efix, emode, alatt, angdeg,...
%                                               u, v, psi, omega, dpsi, gl, gs, grid_size_in, urange_in)
%
% Input: (in the following, nfile = no. spe files)
%   dummy           Dummy sqw object  - used only to ensure that this service routine was called
%   spe_file        Full file name of spe file - character string or cell array of
%                  character strings for more than one file
%   par_file        Full file name of detector parameter file (Tobyfit format)
%   sqw_file        Full file name of output sqw file
%   efix            Fixed energy (meV)                 [scalar or vector length nfile]
%   emode           Direct geometry=1, indirect geometry=2    [scalar]
%   alatt           Lattice parameters (Ang^-1)        [row or column vector]
%   angdeg          Lattice angles (deg)               [row or column vector]
%   u               First vector (1x3) defining scattering plane (r.l.u.)
%   v               Second vector (1x3) defining scattering plane (r.l.u.)
%   psi             Angle of u w.r.t. ki (deg)         [scalar or vector length nfile]
%   psi_planned     As psi, but is a vector specifying the full range of planned measurements
%   omega           Angle of axis of small goniometer arc w.r.t. notional u (deg) [scalar or vector length nfile]
%   dpsi            Correction to psi (deg)            [scalar or vector length nfile]
%   gl              Large goniometer arc angle (deg)   [scalar or vector length nfile]
%   gs              Small goniometer arc angle (deg)   [scalar or vector length nfile]
%   grid_size_in    [Optional] Scalar or row vector of grid dimensions. Default is [50,50,50,50]
%   urange_in       [Optional] Range of data grid for output. If not given, then uses smallest hypercuboid
%                  that encloses the whole data range that would be available if psi_planned were used.
%
% Output:
% --------
%   tmp_file        List of temporary files
%   grid_size       Actual size of grid used (size is unity along dimensions
%                  where there is zero range of the data points)
%   urange          Actual range of grid
%%   Overloaded methods:
%      sqw/accumulate_sqw
%      sqw/accumulate_sqw
%