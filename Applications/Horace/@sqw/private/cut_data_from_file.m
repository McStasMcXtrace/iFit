% Accumulates pixels into bins defined by cut parameters
%
%   >> [s, e, npix, npix_retain] = cut_data (fid, nstart, nend, urange_step, rot_ustep, trans_bott_left, ebin, trans_elo, pax, nbin, keep_pix)
%
% Input:
%   fid             File identifier, with current position in the file being the start of the array of pixel information
%   nstart          Column vector of read start locations in file
%   nend            Column vector of read end locations in file
%   keep_pix        Set to true if wish to retain the information about individual pixels; set to false if not
%   pix_tmpfile_ok  if keep_pix=false, ignore
%                   if keep_pix=true, set buffering option:
%                       pix_tmpfile_ok = false: Require that output argument pix is 9xn array of u1,u2,u3,u4,irun,idet,ien,s,e
%                                              for each retained pixel
%                       pix_tmpfile_ok = true:  Buffering of pixel info to temporary files if pixels exceed a threshold
%                                              In this case, output argument pix contains details of temporary files (see below)
%   urange_step     [2x4] array of the ranges of the data as defined by (i) output proj. axes ranges for
%                  integration axes (or plot axes with one bin), and (ii) step range (0 to no. bins)
%                  for plotaxes (with more than one bin)
%   rot_ustep       Matrix [3x3]     --|  that relate a vector expressed in the
%   trans_bott_left Translation [3x1]--|  frame of the pixel data to no. steps from lower data limit
%                                             r_step(i) = A(i,j)(r(j) - trans(j))
%   ebin            Energy bin width (plays role of rot_ustep for energy axis)
%   trans_elo       Bottom of energy scale (plays role of trans_bott_left for energy axis)
%   pax             Indices of plot axes (with two or more bins) [row vector]
%   nbin            Number of bins along the projection axes with two or more bins [row vector]
%
% Output:
%   s               Array of accumulated signal from all contributing pixels (dimensions match the plot axes)
%   e               Array of accumulated variance
%   npix            Array of number of contributing pixels (if keep_pix==true, otherwise pix=[])
%   urange_step_pix Actual range of contributing pixels
%   pix             if keep_pix=false, pix=[];
%                   if keep_pix==true, then contents depend on value of pix_tmpfile_ok:
%                       pix_tmpfile_ok = false: contains u1,u2,u3,u4,irun,idet,ien,s,e for each retained pixel
%                       pix_tmpfile_ok = true: structure with fields
%                           pix.tmpfiles        cell array of filename(s) containing npix and pix
%                           pix.pos_npixstart   array with position(s) from start of file(s) of array npix
%                           pix.pos_pixstart    array with position(s) from start of file(s) of array npix
%   npix_retain     Number of pixels that contribute to the cut
%   npix_read       Number of pixels read from file
%
%
% Note:
% - Redundant input variables in that urange_step(2,pax)=nbin in implementation of 19 July 2007
% - Aim to take advantage of in-place working within accumulate_cut
%