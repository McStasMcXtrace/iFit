% Create an output sqw file with dummy data using array(s) of energy bins instead spe file(s).
%
%   >> fake_sqw (sqw, en, par_file, sqw_file, efix, emode, alatt, angdeg,...
%                    u, v, psi, omega, dpsi, gl, gs)
%
%   >> fake_sqw (sqw, en, par_file, sqw_file, efix, emode, alatt, angdeg,...
%                    u, v, psi, omega, dpsi, gl, gs, grid_size_in, urange_in)
%
%   >> [tmp_file, grid_size, urange] = fake_sqw (...)
%
% Input:
% ------
%   en              Energy bin boundaries (must be monotonically increasing and equally spaced)
%               or  cell array of arrays of energy bin boundaries, one array per spe file
%   par_file        Full file name of detector parameter file (Tobyfit format)
%   sqw_file        Full file name of output sqw file
%   efix            Fixed energy (meV)                 [scalar or vector length nfile]
%   emode           Direct geometry=1, indirect geometry=2    [scalar]
%   alatt           Lattice parameters (Ang^-1)        [row or column vector]
%   angdeg          Lattice angles (deg)               [row or column vector]
%   u               First vector (1x3) defining scattering plane (r.l.u.)
%   v               Second vector (1x3) defining scattering plane (r.l.u.)
%   psi             Angle of u w.r.t. ki (deg)         [scalar or vector length nfile]
%   omega           Angle of axis of small goniometer arc w.r.t. notional u (deg) [scalar or vector length nfile]
%   dpsi            Correction to psi (deg)            [scalar or vector length nfile]
%   gl              Large goniometer arc angle (deg)   [scalar or vector length nfile]
%   gs              Small goniometer arc angle (deg)   [scalar or vector length nfile]
%   grid_size_in    [Optional] Scalar or row vector of grid dimensions. Default is [50,50,50,50]
%   urange_in       [Optional] Range of data grid for output. If not given, then uses smallest hypercuboid
%                                       that encloses the whole data range.
%
% Output:
% --------
%   tmp_file        List of temporary files
%   grid_size       Actual size of grid used (size is unity along dimensions
%                  where there is zero range of the data points)
%   urange          Actual range of grid
%
%
% Use to generate an sqw file that can be used for creating simulations. Syntax very similar to
% gen_sqw: the only difference is that the input spe data is replaced by energy bin boundaries.
%%   Overloaded methods:
%      sqw/fake_sqw
%      sqw/fake_sqw
%