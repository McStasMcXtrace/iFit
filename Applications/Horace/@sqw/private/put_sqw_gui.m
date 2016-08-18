% Write an sqw data structure to file
%
%   >> [mess, position, npixtot, type] = put_sqw (outfile, main_header, header, detpar, data)
%   >> [mess, position, npixtot, type] = put_sqw (outfile, main_header, header, detpar, data, '-nopix')
%   >> [mess, position, npixtot, type] = put_sqw (outfile, main_header, header, detpar, data, '-pix')
%   >> [mess, position, npixtot, type] = put_sqw (outfile, main_header, header, detpar, data, '-pix',...
%                                                              infiles, npixstart, pixstart, run_label)
%
% Input:
% -------
%   outfile     File name, or file identifier of open file, to which to write data
%   main_header Main header block (for details of data structure, type >> help get_sqw_main_header)
%   header      Header block (for details of data structure, type >> help get_sqw_header)
%   detpar      Detector parameters (for details of data structure, type >> help get_sqw_detpar)
%   data        Valid sqw data structure which must contain the fields listed below  (for details, type >> help get_sqw_data)
%                       type 'b'    fields: uoffset,...,s,e
%                       type 'b+'   fields: uoffset,...,s,e,npix
%                       type 'a'    uoffset,...,s,e,npix,urange,pix
%               In addition, will take the data structure of type 'a' without the individual pixel information ('a-')
%                       type 'a-'   uoffset,...,s,e,npix,urange
%
%   opt         [optional argument for type 'a' or type 'a-'] Determines whether or not to write pixel info, and
%               from which source:
%                 -'-nopix'  do not write the info for individual pixels
%                 -'-pix'    write pixel information
%               The default source of pixel information is the data structure, but if the optional arguments below
%               are given, then use the corresponding source of pixel information
%
%               Can also choose to write just the headaer information in data:
%                 -'-h'      the information as read with '-h' option in get_sqw is written
%                           namely the fields: uoffset,...,dax
%                           (Note: urange will not be written, even if present - types 'a' or 'a-')
%
%   infiles     Cell array of file names, or array of file identifiers of open file, from
%                                   which to accumulate the pixel information
%   npixstart   Position (in bytes) from start of file of the start of the field npix
%   pixstart    Position (in bytes) from start of file of the start of the field pix
%   run_label   Indicates how to re-label the run index (pix(5,...) 
%                       'fileno'    relabel run index as the index of the file in the list infiles
%                       'nochange'  use the run index as in the input file
%                   This option exists to deal with the two limiting cases 
%                    (1) There is one file per run, and the run index in the header block is the file
%                       index e.g. as in the creating of the master sqw file
%                    (2) The run index is already written to the files correctly indexed into the header
%                       e.g. as when temporary files have been written during cut_sqw
%
%
% Output:
% --------
%   mess        If no problem, then mess=''
%               If a problems contains error message and position=[], npixtot=[]; file left open if passed as a fid
%   position    Position (in bytes from start of file) of blocks of fields and large fields:
%                   position.main_header    start of main_header block (=[] if not written)
%                   position.header         start of each header block (header is column vector, length main_header.nfiles)
%                   position.detpar         start of detector parameter block (=[] if not written)
%                   position.data           start of data block
%                   position.s      position of array s
%                   position.e      position of array e
%                   position.npix   position of array npix (=[] if npix not written)
%                   position.pix    position of array pix  (=[] if pix not written)
%   npixtot     Total number of pixels written to file  (=[] if pix not written)
%   type        Type of sqw data written to file: 'a', 'a-', 'b+' or 'b'
% 
%