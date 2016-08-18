% Checks if the detector parameters are equal for two detector structures
% (apart from the names of the files from which they were read)
%
%   >> ans = isequal_par(det1,det2)
%
%   det1, det2      Detector parameter structures of Tobyfit format (see (get_par)
%   ok              Logical, set to true if the structures are equal (apart from the
%                  name and path to the files)
%