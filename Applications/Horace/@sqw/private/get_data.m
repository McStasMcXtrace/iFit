% Load spe file and detector parameter file to create data structure, removing
% masked detector groups and masking any other pixels that contain null data
%
%   >> [data,det,keep,det0]=get_data(spe_file,par_file)
%
%   spe_file        File name of spe file
%   par_file        File name of detector parameter file
%   
%   data            Data structure containing data for unmasked detectors
%                  (see get_spe for list of fields)
%   det             Data structure containing detector parameters for unmasked
%                  detectors (see get_par for fields)
%   keep            List of the detector groups that were unmasked
%   det0            Data structure of detector parameters before masking
%