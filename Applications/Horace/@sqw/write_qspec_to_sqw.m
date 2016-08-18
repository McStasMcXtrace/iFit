% Read ascii column data and create a single sqw file.
%
%   >> write_qspec_to_sqw (dummy, qspec_file, sqw_file, efix, emode, alatt, angdeg,...
%                                                   u, v, psi, omega, dpsi, gl, gs, grid_size_in, urange_in)
%
% Input:
%   dummy           Dummy sqw object  - used only to ensure that this service routine was called
%   qspec_file      Full file name of ascii file containing qx-qy-qz-eps-signal-error column data.
%                   Here qz  is the component of momentum along ki (Ang^-1)
%                        qy  is component vertically upwards (Ang^-1)
%                        qx  defines a hight-hand coordinate frame with qy' and qz'
%                        S   signal
%                        ERR standard deviation
%                       
%   sqw_file        Full file name of output sqw file
%
%   efix            Fixed energy (meV) (if elastic data ie. emode=0, the value will be ignored)
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
%                  that encloses the whole data range
%
% Output:
%   grid_size       Actual grid size used (size is unity along dimensions
%                  where there is zero range of the data points)
%   urange          Actual range of grid
%%   Overloaded methods:
%      sqw/write_qspec_to_sqw
%      sqw/write_qspec_to_sqw
%