% Determine the un-masked detectors in an SPE data file, and a logical array of
% any other pixels to mask in other detectors
%
%   >> [keep,pixmask]=mask_spe(data)
%   
%   data        spe data structure as read by get_spe
%   
%   keep        Row vector of the detector groups that are unmasked
%              These are the groups where all pixels are masked
%   pixmask 	Logical array, size(pixmask)=data.S(:,keep), with true
%              at those pixels with need to masked. This is returned as
%              empty if there are no pixels to be masked.
%