% Load data from ASCII Tobyfit .par file
%   >> par = get_par_matlab(filename)
%
%     filename      name of par file
%
%     par(5,ndet)   contents of array
%
%     1st column    sample-detector distance
%     2nd  "        scattering angle (deg)
%     3rd  "        azimuthal angle (deg)
%                   (west bank = 0 deg, north bank = -90 deg etc.)
%                   (Note the reversed sign convention cf .phx files)
%     4th  "        width (m)
%     5th  "        height (m)
%