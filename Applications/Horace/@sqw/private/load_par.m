% Load data from ASCII Tobyfit .par file
%   >> det = get_par(filename)
%
% data has following fields:
%   det.filename    Name of file excluding path
%   det.filepath    Path to file including terminating file separator
%   det.group       Row vector of detector group number - assumed to be 1:ndet
%   det.x2          Secondary flightpath (m)
%   det.phi         Row vector of scattering angles (deg)
%   det.azim        Row vector of azimuthal angles (deg)
%                  (West bank=0 deg, North bank=90 deg etc.)
%   det.width       Row vector of detector widths (m)
%   det.height      Row vector of detector heights (m)
%
%  if varargin present, do not convert into detector structure but return
%  initial array 
%%   Overloaded methods:
%      loader_nxspe/load_par
%      loader_ascii/load_par
%