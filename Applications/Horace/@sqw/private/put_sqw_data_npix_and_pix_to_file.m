% Write npix and pix to a file with same format as write_sqw_data
%
%   >> [mess, position] = put_sqw_data_npix_and_pix_to_file (outfile, npix, pix)
%
% Input:
% ------
%   outfile     File name, or file identifier of open file, to which to append data
%   npix        Array containing the number of pixels in each bin
%   data.pix    Array containing data for each pixel:
%              If npixtot=sum(npix), then pix(9,npixtot) contains:
%                   u1      -|
%                   u2       |  Coordinates of pixel in the projection axes of the original sqw file(s)
%                   u3       |
%                   u4      -|
%                   irun        Run index in the header block from which pixel came
%                   idet        Detector group number in the detector listing for the pixel
%                   ien         Energy bin number for the pixel in the array in the (irun)th header
%                   signal      Signal array
%                   err         Error array (variance i.e. error bar squared)
%
% Output:
% -------
%   mess        Message if there was a problem writing; otherwise mess=''
%   position    Position (in bytes from start of file) of large fields:
%                   position.npix   position of array npix (in bytes) from beginning of file
%                   position.pix    position of array pix (in bytes) from beginning of file
%   npixtot     Total number of pixels written to file  (=[] if pix not written)
%