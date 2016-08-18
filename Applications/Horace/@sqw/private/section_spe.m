% Obtain reduced spe data file consisting of those bins entirely contained within
% given detector number range and energy transfer range
%
%   >> data_new=section_spe(data,det_array)
%   >> data_new=section_spe(data,det_array,en_range)
%
% Input:
%   data        spe data structure (see get_spe)
%   det_array   [det1,det2,...detn] are detector numbers to keep
%              (Note: the order is retained regardless if monotonic or not)
%   en_range    (Optional): [en_lo, en_hi] is energy transfer range to keep
%              Default is to keep full energy range
%
% Output:
%   data_new    Output spe data structure
%
% e.g.
%   >> dnew=section_spe(d,[15:50],[12,15]);
%