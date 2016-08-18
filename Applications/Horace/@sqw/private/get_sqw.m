% Load an sqw file from disk
%
% Syntax:
%   >> [main_header,header,detpar,data,mess,position,npixtot,type,current_format] = get_sqw (infile)
%   >> [main_header,header,detpar,data,mess,position,npixtot,type,current_format] = get_sqw (infile, '-h')
%   >> [main_header,header,detpar,data,mess,position,npixtot,type,current_format] = get_sqw (infile, '-hverbatim')
%   >> [main_header,header,detpar,data,mess,position,npixtot,type,current_format] = get_sqw (infile, '-nopix')
%   >> [main_header,header,detpar,data,mess,position,npixtot,type,current_format] = get_sqw (infile, npix_lo, npix_hi)
%
% Input:
% --------
%   infile      File name, or file identifier of open file, from which to read data
%   opt         [optional] Determines which fields to read from the data block:
%                   '-h'            Header-type information only: fields read: 
%                                       filename, filepath, title, alatt, angdeg,...
%                                           uoffset,u_to_rlu,ulen,ulabel,iax,iint,pax,p,dax[,urange]
%                                  (If file was written from a structure of type 'b' or 'b+', then
%                                  urange does not exist, and the output field will not be created)
%                   '-hverbatim'    Same as '-h' except that the file name as stored in the main_header and
%                                  data sections are returned as stored, not constructed from the
%                                  value of fopen(fid). This is needed in some applications where
%                                  data is written back to the file with a few altered fields.
%                   '-nopix'        Pixel information not read (only meaningful for sqw data type 'a')
%
%                    Default: read all fields of the corresponding sqw data type ('b','b+','a','a-')
%
%   npix_lo     -|- [optional] pixel number range to be read from the file
%   npix_hi     -|
%
% Output:
% --------
%   main_header Main header block (for details of data structure, type >> help get_sqw_main_header)
%   header      Header block (for details of data structure, type >> help get_sqw_header)
%              Cell array if more than one contributing spe file.
%   detpar      Detector parameters (for details of data structure, type >> help get_sqw_detpar)
%   data        Output data structure which will contain the fields listed below (for details, type >> help get_sqw_data) 
%                       type 'b'    fields: filename,...,dax,s,e
%                       type 'b+'   fields: filename,...,dax,s,e,npix
%                       type 'a'    fields: filename,...,dax,s,e,npix,urange,pix
%                       type 'a-'   fields: filename,...,dax,s,e,npix,urange
%               or header information
%   mess        Error message; blank if no errors, non-blank otherwise
%   position    Position (in bytes from start of file) of blocks of fields and large fields:
%                   position.main_header    start of main_header block (=[] if not written)
%                   position.header         start of each header block (header is column vector, length main_header.nfiles)
%                   position.detpar         start of detector parameter block (=[] if not written)
%                   position.data           start of data block
%                   position.s      position of array s
%                   position.e      position of array e
%                   position.npix   position of array npix (=[] if npix not written)
%                   position.pix    position of array pix  (=[] if pix not written)
%   npixtot     Total number of pixels written to file (=[] if pix not present)
%   type        Type of sqw data written to file: 
%               Valid sqw data structure, which must contain the fields listed below 
%                       type 'b'    fields: filename,...,dax,s,e
%                       type 'b+'   fields: filename,...,dax,s,e,npix
%                       type 'a'    fields: filename,...,dax,s,e,npix,urange,pix
%               or if the pix field is not read from type 'a', in which case 
%                       type 'a-'   fields: filename,...,dax,s,e,npix,urange
%