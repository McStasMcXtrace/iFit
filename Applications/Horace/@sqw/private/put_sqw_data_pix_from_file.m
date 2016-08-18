% Write pixel information to file, reading that pixel information from a collection of other files
%
% Syntax:
%   >> mess = put_sqw_data_pix_from_file (fid, infiles, npixstart, pixstart)
%
% Input:
%   fout            File identifier of output file (opened for binary writing)
%   infiles         Cell array of file names, or array of file identifiers of open files, from
%                  which to accumulate the pixel information
%   pos_npixstart   Position (in bytes) from start of file of the start of the field npix
%   pos_pixstart    Position (in bytes) from start of file of the start of the field pix
%   npix_cumsum     Accumulated sum of number of pixels per bin across all the files
%   run_label       Indicates how to re-label the run index (pix(5,...) 
%                       'fileno'    relabel run index as the index of the file in the list infiles
%                       'nochange'  use the run index as in the input file
%                   This option exists to deal with the two limiting cases 
%                    (1) There is one file per run, and the run index in the header block is the file
%                       index e.g. as in the creating of the master sqw file
%                    (2) The run index is already written to the files correctly indexed into the header
%                       e.g. as when temporary files have been written during cut_sqw
%
% Output:
%   mess            Message if there was a problem writing; otherwise mess=''
%
% Notes:
%   Take care when using this function. No checks are performed that the input files have the
%  correct length of arrays npix and pix. It is assumed that this checking has already been done.
%
%  The reason for this function is that the output sqw structure may be too large to be held in memory.
% This happens in particular during construction of the 'master' sqw file from a collection of sqw files, and
% from taking large cuts from an sqw file (during which temporary files are written with the pixel information to
% avoid out-of-memory problems).
% 
%