% Read a single spe file and a detector parameter file, and create a single sqw file.
% to file.
%
%   >> write_spe_to_sqw (dummy, spe_file, par_file, sqw_file, efix, emode, alatt, angdeg,...
%                                                   u, v, psi, omega, dpsi, gl, gs, grid_size_in, urange_in)
%
% Input:
%   dummy           Dummy sqw object  - used only to ensure that this service routine was called
%   spe_data        Source of spe data e.g. full file name of spe file or nxspe file
%   par_file        Full file name of detector parameter file (Tobyfit format)
%   sqw_file        Full file name of output sqw file
%
%   efix            Fixed energy (meV) (if elastic data ie. emode=0, the value will be ignored and set to zero internally)
%   emode           Direct geometry=1, indirect geometry=2, elastic=0
%   alatt           Lattice parameters (Ang^-1)
%   angdeg          Lattice angles (deg)
%   u               First vector (1x3) defining scattering plane (r.l.u.)
%   v               Second vector (1x3) defining scattering plane (r.l.u.)
%   psi             Angle of u w.r.t. ki (rad)
%   omega           Angle of axis of small goniometer arc w.r.t. notional u
%   dpsi            Correction to psi (rad)
%   gl              Large goniometer arc angle (rad)
%   gs              Small goniometer arc angle (rad)
%   grid_size_in    Scalar or row vector of grid dimensions. Default is [1x1x1x1]
%   urange_in       Range of data grid for output. If not given, then uses smallest hypercuboid
%                   that encloses the whole data range
%   det0            if present, refers to the array of unmasked detector parameters in horace data format
%                   the array is written in sqw file as reference detectors
%                   allowing to combine multiple sqw (or tmp) files togeter
%                   as masked detectors are usually different for different
%                   runs. If it is absent, appropriate data file is taken
%                   from the par_file or rundata
%
% Output:
%   grid_size       Actual grid size used (size is unity along dimensions
%                  where there is zero range of the data points)
%   urange          Actual range of grid
%%   Overloaded methods:
%      sqw/write_spe_to_sqw
%      sqw/write_spe_to_sqw
%