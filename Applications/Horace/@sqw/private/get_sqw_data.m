% Read the data block from an sqw file.
% The file pointer is left at the end of the data block.
%
% Syntax:
%   >> [data, mess] = get_sqw_data(fid)
%   >> [data, mess] = get_sqw_data(fid, opt)
%   >> [data, mess] = get_sqw_data(fid, npix_lo, npix_hi)
%
% To any of the above: to read original prototype format file
%   >> [data, mess] = get_sqw_data(..., '-prototype')
%
%   * The fields filename, filepath, title, alatt, and angdeg are not stored in this prototype format,
%     Fields filename and filepath are constructed in this routine, but title, alatt and angdeg are
%     just given dummy values which must be filled from the main_header and header information.
%
%   * This option is not compatible with files of type 'b', because the npix information is needed to convert
%     the signal and error in each bin into the new format by normalising by number of pixels in each bin.
%
% Input:
% ------
%   fid         File pointer to (already open) binary file
%   data_in     [optional] Data structure to which the data fields below will be added or overwrite.
%   opt         [optional] Determines which fields to read
%                   '-h'     header-type information only: fields read: 
%                               filename, filepath, title, alatt, angdeg,...
%                                   uoffset,u_to_rlu,ulen,ulabel,iax,iint,pax,p,dax[,urange]
%                              (If file was written from a structure of type 'b' or 'b+', then
%                               urange does not exist, and the output field will not be created)
%                   '-hverbatim'    Same as '-h' except that the file name as stored in the main_header and
%                                  data sections are returned as stored, not constructed from the
%                                  value of fopen(fid). This is needed in some applications where
%                                  data is written back to the file with a few altered fields.
%                   '-nopix' Pixel information not read (only meaningful for sqw data type 'a')
%
%                    Default: read all fields of the corresponding sqw data type ('b','b+','a','a-')
%
%   npix_lo     -|- [optional] pixel number range to be read from the file 
%   npix_hi     -|
%
% Output:
% -------
%   data        Output data structure which must contain the fields listed below 
%                       type 'b'    fields: filename,...,dax,s,e
%                       type 'b+'   fields: filename,...,dax,s,e,npix
%                       type 'a'    fields: filename,...,dax,s,e,npix,urange,pix
%                       type 'a-'   fields: filename,...,dax,s,e,npix,urange
%               or header information   
%   mess        Error message; blank if no errors, non-blank otherwise
%   position    Position (in bytes from start of file) of large fields:
%                   position.s      position of array s
%                   position.e      position of array e
%                   position.npix   position of array npix (=[] if npix not present)
%                   position.pix    position of array pix (=[] if pix not present)
%   npixtot     Total number of pixels written to file (=[] if pix not present)
%   type        Type of sqw data written to file: 
%               Valid sqw data structure, which will contain the fields listed below 
%                       type 'b'    fields: filename,...,dax,s,e
%                       type 'b+'   fields: filename,...,dax,s,e,npix
%                       type 'a'    fields: filename,...,dax,s,e,npix,urange,pix
%               or if the pix field is not read from type 'a', in which case 
%                       type 'a-'   fields: filename,...,dax,s,e,npix,urange
%
%
% Fields read from the file are:
%
%   data.filename   Name of sqw file that is being read, excluding path
%   data.filepath   Path to sqw file that is being read, including terminating file separator
%          [Note that the filename and filepath that are written to file are ignored; we fill with the 
%           values corresponding to the file that is being read.]
%
%   data.title      Title of sqw data structure
%   data.alatt      Lattice parameters for data field (Ang^-1)
%   data.angdeg     Lattice angles for data field (degrees)
%   data.uoffset    Offset of origin of projection axes in r.l.u. and energy ie. [h; k; l; en] [column vector]
%   data.u_to_rlu   Matrix (4x4) of projection axes in hkle representation
%                      u(:,1) first vector - u(1:3,1) r.l.u., u(4,1) energy etc.
%   data.ulen       Length of projection axes vectors in Ang^-1 or meV [row vector]
%   data.ulabel     Labels of the projection axes [1x4 cell array of character strings]
%   data.iax        Index of integration axes into the projection axes  [row vector]
%                  Always in increasing numerical order
%                       e.g. if data is 2D, data.iax=[1,3] means summation has been performed along u1 and u3 axes
%   data.iint       Integration range along each of the integration axes. [iint(2,length(iax))]
%                       e.g. in 2D case above, is the matrix vector [u1_lo, u3_lo; u1_hi, u3_hi]
%   data.pax        Index of plot axes into the projection axes  [row vector]
%                  Always in increasing numerical order
%                       e.g. if data is 3D, data.pax=[1,2,4] means u1, u2, u4 axes are x,y,z in any plotting
%                                       2D, data.pax=[2,4]     "   u2, u4,    axes are x,y   in any plotting
%   data.p          Cell array containing bin boundaries along the plot axes [column vectors]
%                       i.e. row cell array{data.p{1}, data.p{2} ...} (for as many plot axes as given by length of data.pax)
%   data.dax        Index into data.pax of the axes for display purposes. For example we may have 
%                  data.pax=[1,3,4] and data.dax=[3,1,2] This means that the first plot axis is data.pax(3)=4,
%                  the second is data.pax(1)=1, the third is data.pax(2)=3. The reason for data.dax is to allow
%                  the display axes to be permuted but without the contents of the fields p, s,..pix needing to
%                  be reordered [row vector]
%   data.s          Cumulative signal.  [size(data.s)=(length(data.p1)-1, length(data.p2)-1, ...)]
%   data.e          Cumulative variance [size(data.e)=(length(data.p1)-1, length(data.p2)-1, ...)]
%   data.npix       No. contributing pixels to each bin of the plot axes.
%                  [size(data.pix)=(length(data.p1)-1, length(data.p2)-1, ...)]
%   data.urange     True range of the data along each axis [urange(2,4)]
%   data.pix        Array containing data for eaxh pixel:
%                  If npixtot=sum(npix), then pix(9,npixtot) contains:
%                   u1      -|
%                   u2       |  Coordinates of pixel in the projection axes
%                   u3       |
%                   u4      -|
%                   irun        Run index in the header block from which pixel came
%                   idet        Detector group number in the detector listing for the pixel
%                   ien         Energy bin number for the pixel in the array in the (irun)th header
%                   signal      Signal array
%                   err         Error array (variance i.e. error bar squared)
%
%
% Notes:
% ------
%   It is assumed that the file corresponds to a valid type (i.e. that any use with implementation of sqw as
%   a proper object has already checked the consistency of the fields).
%