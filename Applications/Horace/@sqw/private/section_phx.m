% Obtain reduced phx data file for a given subsection
%
%   >> data_new=section_phx(data,det_array)
%
% Input:
%   det         phx data structure (see get_phx)
%   det_array   [det1,det2,...detn] is a list of the indexes of the entries in
%              detector arrays that are to be kept. THis is only the same as the 
%              detector group numbers if det.group = 1:ndet
%              (Note: the order is retained regardless if monotonic or not)
% Ouput:
%   det_new    Output phx data structure
%