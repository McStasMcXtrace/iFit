function data = read_nii(filename)
% read a NifTi Neuroimaging Informatics Technology Initiative volume
%  data = read_nii(filename)
%
% <http://nifti.nimh.nih.gov/>
%
% references:
%   http://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
%   Jimmy Shen 23 Oct 2005 (Updated 22 Jan 2014)
% and
%   <http://www.mathworks.com/matlabcentral/fileexchange/29344-read-medical-data-3d>
%   by Dirk-Jan Kroon 10 Nov 2010 (Updated 23 Feb 2011)
%
% See also: read_poscar, read_pdb, read_xyz, read_mrc, read_analyze

data = [];
if nargin == 0, return; end

try
  % a good, reliable loader by Jimmy Shen
  data = load_nii(filename);
catch
  % an older loader by Dirk-Jan Kroon
  info = nii_read_header(filename);
  if ~isempty(info)
    data.Header = info;
    data.Image  = nii_read_volume(info);
  end
end

% ------------------------------------------------------------------------------
%Copyright (c) 2010, Dirk-Jan Kroon
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are
%met:
%
%* Redistributions of source code must retain the above copyright
%  notice, this list of conditions and the following disclaimer.
%* Redistributions in binary form must reproduce the above copyright
%  notice, this list of conditions and the following disclaimer in
%  the documentation and/or other materials provided with the distribution
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%POSSIBILITY OF SUCH DAMAGE.

function [info] = nii_read_header(filename)
% function for reading header of NifTi ( .nii ) volume file
%
% info = nii_read_header(filename);
%
% examples:
% 1,  info=nii_read_header()
% 2,  info=nii_read_header('volume.nii');

if(exist('filename','var')==0)
     [filename, pathname] = uigetfile('*.nii', 'Read nii-file');
     filename = [pathname filename];
end

bswap=false;
test=true;
while(test)
    if(bswap)
        fid=fopen(filename,'rb','b');
    else
        fid=fopen(filename,'rb','l');
    end    
    if(fid<0)
         fprintf('could not open file %s\n',filename);
         return
    end

    %get the file size
    fseek(fid,0,'eof');
    info.Filesize = ftell(fid); 
    fseek(fid,0,'bof');
    info.Filename=filename;
    info.SizeofHdr=fread(fid,1,'int');
    info.DataType=fread(fid, 10, 'uint8=>char')';
    info.DbName=fread(fid, 18, 'uint8=>char')';
    info.Extents=fread(fid,1,'int');
    info.SessionError=fread(fid,1,'uint16');
    info.Regular=fread(fid, 1, 'uint8=>char')';
    info.DimInfo=fread(fid, 1, 'uint8=>char')';
    swaptemp=fread(fid, 1, 'uint16')';
    info.Dimensions=fread(fid,7,'uint16')'; % dim = [ number of dimensions x,y,z,t,c1,c2,c3];
	
    if(swaptemp(1)<1||swaptemp(1)>7), bswap=true; fclose(fid); else test=false; end
    
    if info.SizeofHdr > info.Filesize || info.Extents > info.Filesize
      fprintf([ mfilename ': wrong file hdr extent %i %i\n' ], info.SizeofHdr, info.Extents);
      info = [];
      fclose(fid);
      return
    end
end
info.headerbswap=bswap;
info.IntentP1=fread(fid,1,'float');
info.IntentP2=fread(fid,1,'float');
info.IntentP3=fread(fid,1,'float');
info.IntentCode=fread(fid,1,'uint16');
info.DataType=fread(fid,1,'uint16');
datatypestr{1}={0,'UNKNOWN',  0}; % what it says, dude           
datatypestr{2}={1,'BINARY',   1}; % binary (1 bit/voxel)         
datatypestr{3}={2,'UINT8'  ,  8};% unsigned char (8 bits/voxel) 
datatypestr{4}={4,'INT16'   , 16}; % signed short (16 bits/voxel) 
datatypestr{5}={8,'INT32'  ,  32}; % signed int (32 bits/voxel)   
datatypestr{6}={16,'FLOAT' ,  32}; % float (32 bits/voxel)        
datatypestr{7}={32,'COMPLEX', 64}; % complex (64 bits/voxel)      
datatypestr{8}={64,'DOUBLE',  64}; % double (64 bits/voxel)       
datatypestr{9}={128,'RGB'  ,  24}; % RGB triple (24 bits/voxel)   
datatypestr{10}={255,'ALL'  ,  0}; % not very useful (?)          
datatypestr{11}={256,'INT8' ,  8}; % signed char (8 bits)         
datatypestr{12}={512,'UINT16', 16}; % unsigned short (16 bits)     
datatypestr{13}={768,'UINT32', 32}; % unsigned int (32 bits)       
datatypestr{14}={1024,'INT64', 64}; % long long (64 bits)          
datatypestr{15}={1280,'UINT64',     64}; % unsigned long long (64 bits) 
datatypestr{16}={1536,'FLOAT128',   128}; % long double (128 bits)       
datatypestr{17}={1792,'COMPLEX128', 128}; % double pair (128 bits)       
datatypestr{18}={2048,'COMPLEX256', 256}; % long double pair (256 bits)  
datatypestr{19}={2304,'RGBA32', 32}; % 4 byte RGBA (32 bits/voxel) 
info.datatypestr='UNKNOWN';
info.bitvoxel=0;
for i=1:19
    if(datatypestr{i}{1}==info.DataType)
        info.DataTypeStr=datatypestr{i}{2};
        info.BitVoxel=datatypestr{i}{3};
    end
end

info.Bitpix=fread(fid,1,'uint16');
info.SliceStart=fread(fid,1,'uint16');
temp=fread(fid,1,'float');
info.PixelDimensions=fread(fid,7,'float');

info.VoxOffset=fread(fid,1,'float');
info.RescaleSlope=fread(fid,1,'float');
info.RescaleIntercept=fread(fid,1,'float');
info.SliceEnd=fread(fid,1,'uint16');
info.SliceCode=fread(fid, 1, 'uint8=>char')';
info.XyztUnits=fread(fid, 1, 'uint8')';
dataunitsstr{1}={'UNKNOWN', 0}; %! NIFTI code for unspecified units. 
dataunitsstr{2}={'METER',   1};  %! NIFTI code for meters. 
dataunitsstr{3}={'MM',    2};  %! NIFTI code for millimeters. 
dataunitsstr{4}={'MICRON ', 3};  %! NIFTI code for micrometers. 
dataunitsstr{5}={'SEC',    8};  %! NIFTI code for seconds. 
dataunitsstr{6}={'MSEC',   16};  %! NIFTI code for milliseconds. 
dataunitsstr{7}={'USEC',  24};  %! NIFTI code for microseconds. 
dataunitsstr{8}={'HZ',  32};  %! NIFTI code for Hertz. 
dataunitsstr{9}={'PPM',  40};  %! NIFTI code for ppm. 
dataunitsstr{10}={'RADS',  48};  %! NIFTI code for radians per second. 
info.xyzt_unitsstr='UNKNOWN';
for i=1:10,
    if(dataunitsstr{i}{2}==info.XyztUnits)
        info.XyztUnitsStr=dataunitsstr{i}{1};
    end
end

info.CalMax=fread(fid,1,'float');
info.CalMin=fread(fid,1,'float');
info.Slice_duration=fread(fid,1,'float');
info.Toffset=fread(fid,1,'float');
info.Glmax=fread(fid,1,'int');
info.Glmin=fread(fid,1,'int');
info.Descrip=fread(fid, 80, 'uint8=>char')';
info.AuxFile=fread(fid, 24, 'uint8=>char')';
info.QformCode=fread(fid,1,'uint16');
info.SformCode=fread(fid,1,'uint16');
info.QuaternB=fread(fid,1,'float');
info.QuaternC=fread(fid,1,'float');
info.QuaternD=fread(fid,1,'float');
info.QoffsetX=fread(fid,1,'float');
info.QoffsetY=fread(fid,1,'float');
info.QoffsetZ=fread(fid,1,'float');
info.SrowX=fread(fid,4,'float');
info.SrowY=fread(fid,4,'float');
info.SrowZ=fread(fid,4,'float');
info.IntentName=fread(fid, 16, 'uint8=>char')';
info.Magic=fread(fid, 4, 'uint8=>char')';

fclose(fid);
if ~strncmp(info.Magic, 'nii', 3) && ~strncmp(info.Magic, 'n+1', 3)
  error([ mfilename ': ' filename ' is probably not a NifTI volume file.' ])
end

% ------------------------------------------------------------------------------
function V = nii_read_volume(info)
% function for reading volume of NifTi ( .nii ) volume file
% nii_read_header(file-info)
%
% volume = nii_read_volume(file-header)
%
% examples:
% 1: info = nii_read_header()
%    V = nii_read_volume(info);
%    imshow(squeeze(V(:,:,round(end/2))),[]);
%
% 2: V = nii_read_volume('test.nii');

if(~isstruct(info)) info=nii_read_header(info); end

% Open v3d file
fid=fopen(info.Filename,'rb');

  % Seek volume data start
  datasize=prod(info.Dimensions)*(info.BitVoxel/8);
  fsize=info.Filesize;
  fseek(fid,fsize-datasize,'bof');

  % Read Volume data
  switch(info.DataTypeStr)
      case 'INT8'
        V = int8(fread(fid,datasize,'int8'));
      case 'UINT8'
        V = uint8(fread(fid,datasize,'uint8'));
      case 'INT16'
        V = int16(fread(fid,datasize,'int16'));
      case 'UINT16'
        V = uint16(fread(fid,datasize,'uint16'));
      case 'INT32'
        V = int32(fread(fid,datasize,'int32'));
      case 'UINT32'
        V = uint32(fread(fid,datasize,'uint32'));
      case 'INT64'
        V = int64(fread(fid,datasize,'int64'));
      case 'UINT64'
        V = uint64(fread(fid,datasize,'uint64'));    
      otherwise
        V = uint8(fread(fid,datasize,'uint8'));
  end
fclose(fid);

% Reshape the volume data to the right dimensions
V = reshape(V,info.Dimensions);





















% ------------------------------------------------------------------------------
%Copyright (c) 2014, Jimmy Shen
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are
%met:
%
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in
%      the documentation and/or other materials provided with the distribution
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%POSSIBILITY OF SUCH DAMAGE.
%
% http://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
% Jimmy Shen 23 Oct 2005 (Updated 22 Jan 2014)


%  Load NIFTI or ANALYZE dataset. Support both *.nii and *.hdr/*.img
%  file extension. If file extension is not provided, *.hdr/*.img will
%  be used as default.
%
%  A subset of NIFTI transform is included. For non-orthogonal rotation,
%  shearing etc., please use 'reslice_nii.m' to reslice the NIFTI file.
%  It will not cause negative effect, as long as you remember not to do
%  slice time correction after reslicing the NIFTI file. Output variable
%  nii will be in RAS orientation, i.e. X axis from Left to Right,
%  Y axis from Posterior to Anterior, and Z axis from Inferior to
%  Superior.
%  
%  Usage: nii = load_nii(filename, [img_idx], [dim5_idx], [dim6_idx], ...
%			[dim7_idx], [old_RGB], [tolerance], [preferredForm])
%  
%  filename  - 	NIFTI or ANALYZE file name.
%  
%  img_idx (optional)  -  a numerical array of 4th dimension indices,
%	which is the indices of image scan volume. The number of images
%	scan volumes can be obtained from get_nii_frame.m, or simply
%	hdr.dime.dim(5). Only the specified volumes will be loaded. 
%	All available image volumes will be loaded, if it is default or
%	empty.
%
%  dim5_idx (optional)  -  a numerical array of 5th dimension indices.
%	Only the specified range will be loaded. All available range
%	will be loaded, if it is default or empty.
%
%  dim6_idx (optional)  -  a numerical array of 6th dimension indices.
%	Only the specified range will be loaded. All available range
%	will be loaded, if it is default or empty.
%
%  dim7_idx (optional)  -  a numerical array of 7th dimension indices.
%	Only the specified range will be loaded. All available range
%	will be loaded, if it is default or empty.
%
%  old_RGB (optional)  -  a scale number to tell difference of new RGB24
%	from old RGB24. New RGB24 uses RGB triple sequentially for each
%	voxel, like [R1 G1 B1 R2 G2 B2 ...]. Analyze 6.0 from AnalyzeDirect
%	uses old RGB24, in a way like [R1 R2 ... G1 G2 ... B1 B2 ...] for
%	each slices. If the image that you view is garbled, try to set 
%	old_RGB variable to 1 and try again, because it could be in
%	old RGB24. It will be set to 0, if it is default or empty.
%
%  tolerance (optional) - distortion allowed in the loaded image for any
%	non-orthogonal rotation or shearing of NIfTI affine matrix. If 
%	you set 'tolerance' to 0, it means that you do not allow any 
%	distortion. If you set 'tolerance' to 1, it means that you do 
%	not care any distortion. The image will fail to be loaded if it
%	can not be tolerated. The tolerance will be set to 0.1 (10%), if
%	it is default or empty.
%
%  preferredForm (optional)  -  selects which transformation from voxels
%	to RAS coordinates; values are s,q,S,Q.  Lower case s,q indicate
%	"prefer sform or qform, but use others if preferred not present". 
%	Upper case indicate the program is forced to use the specificied
%	tranform or fail loading.  'preferredForm' will be 's', if it is
%	default or empty.	- Jeff Gunter
%
%  Returned values:
%  
%  nii structure:
%
%	hdr -		struct with NIFTI header fields.
%
%	filetype -	Analyze format .hdr/.img (0); 
%			NIFTI .hdr/.img (1);
%			NIFTI .nii (2)
%
%	fileprefix - 	NIFTI filename without extension.
%
%	machine - 	machine string variable.
%
%	img - 		3D (or 4D) matrix of NIFTI data.
%
%	original -	the original header before any affine transform.
%  
%  Part of this file is copied and modified from:
%  http://www.mathworks.com/matlabcentral/fileexchange/1878-mri-analyze-tools
%  
%  NIFTI data format can be found on: http://nifti.nimh.nih.gov
%  
%  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)
%
function nii = load_nii(filename, img_idx, dim5_idx, dim6_idx, dim7_idx, ...
			old_RGB, tolerance, preferredForm)

   if ~exist('filename','var')
      error('Usage: nii = load_nii(filename, [img_idx], [dim5_idx], [dim6_idx], [dim7_idx], [old_RGB], [tolerance], [preferredForm])');
   end

   if ~exist('img_idx','var') | isempty(img_idx)
      img_idx = [];
   end

   if ~exist('dim5_idx','var') | isempty(dim5_idx)
      dim5_idx = [];
   end

   if ~exist('dim6_idx','var') | isempty(dim6_idx)
      dim6_idx = [];
   end

   if ~exist('dim7_idx','var') | isempty(dim7_idx)
      dim7_idx = [];
   end

   if ~exist('old_RGB','var') | isempty(old_RGB)
      old_RGB = 0;
   end

   if ~exist('tolerance','var') | isempty(tolerance)
      tolerance = 0.1;			% 10 percent
   end

   if ~exist('preferredForm','var') | isempty(preferredForm)
      preferredForm= 's';		% Jeff
   end

   %  Read the dataset header
   %
   [nii.hdr,nii.filetype,nii.fileprefix,nii.machine] = load_nii_hdr(filename);

   %  Read the header extension
   %
%   nii.ext = load_nii_ext(filename);

   %  Read the dataset body
   %
   [nii.img,nii.hdr] = load_nii_img(nii.hdr,nii.filetype,nii.fileprefix, ...
		nii.machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB);

   %  Perform some of sform/qform transform
   %
   nii = xform_nii(nii, tolerance, preferredForm);

   return					% load_nii

%  internal function

%  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)

function [hdr, filetype, fileprefix, machine] = load_nii_hdr(fileprefix)

   if ~exist('fileprefix','var'),
      error('Usage: [hdr, filetype, fileprefix, machine] = load_nii_hdr(filename)');
   end

   machine = 'ieee-le';
   new_ext = 0;

   if findstr('.nii',fileprefix)
      new_ext = 1;
      fileprefix = strrep(fileprefix,'.nii','');
   end

   if findstr('.hdr',fileprefix)
      fileprefix = strrep(fileprefix,'.hdr','');
   end

   if findstr('.img',fileprefix)
      fileprefix = strrep(fileprefix,'.img','');
   end

   if new_ext
      fn = sprintf('%s.nii',fileprefix);

      if ~exist(fn)
         msg = sprintf('Cannot find file "%s.nii".', fileprefix);
         error(msg);
      end
   else
      fn = sprintf('%s.hdr',fileprefix);

      if ~exist(fn)
         msg = sprintf('Cannot find file "%s.hdr".', fileprefix);
         error(msg);
      end
   end

   fid = fopen(fn,'r',machine);
    
   if fid < 0,
      msg = sprintf('Cannot open file %s.',fn);
      error(msg);
   else
      fseek(fid,0,'bof');

      if fread(fid,1,'int32') == 348
         hdr = read_header(fid);
         fclose(fid);
      else
         fclose(fid);

         %  first try reading the opposite endian to 'machine'
         %
         switch machine,
         case 'ieee-le', machine = 'ieee-be';
         case 'ieee-be', machine = 'ieee-le';
         end

         fid = fopen(fn,'r',machine);

         if fid < 0,
            msg = sprintf('Cannot open file %s.',fn);
            error(msg);
         else
            fseek(fid,0,'bof');

            if fread(fid,1,'int32') ~= 348

               %  Now throw an error
               %
               msg = sprintf('File "%s" is corrupted.',fn);
               fclose(fid);
               error(msg);
            end

            hdr = read_header(fid);
            fclose(fid);
         end
      end
   end

   if strcmp(hdr.hist.magic, 'n+1')
      filetype = 2;
   elseif strcmp(hdr.hist.magic, 'ni1')
      filetype = 1;
   else
      filetype = 0;
   end

   return					% load_nii_hdr


%---------------------------------------------------------------------
function [ dsr ] = read_header(fid)

        %  Original header structures
	%  struct dsr
	%       { 
	%       struct header_key hk;            /*   0 +  40       */
	%       struct image_dimension dime;     /*  40 + 108       */
	%       struct data_history hist;        /* 148 + 200       */
	%       };                               /* total= 348 bytes*/

    dsr.hk   = header_key(fid);
    dsr.dime = image_dimension(fid);
    dsr.hist = data_history(fid);

    %  For Analyze data format
    %
    if ~strcmp(dsr.hist.magic, 'n+1') & ~strcmp(dsr.hist.magic, 'ni1')
        dsr.hist.qform_code = 0;
        dsr.hist.sform_code = 0;
    end

    return					% read_header


%---------------------------------------------------------------------
function [ hk ] = header_key(fid)

    fseek(fid,0,'bof');
    
	%  Original header structures	
	%  struct header_key                     /* header key      */ 
	%       {                                /* off + size      */
	%       int sizeof_hdr                   /*  0 +  4         */
	%       char data_type[10];              /*  4 + 10         */
	%       char db_name[18];                /* 14 + 18         */
	%       int extents;                     /* 32 +  4         */
	%       short int session_error;         /* 36 +  2         */
	%       char regular;                    /* 38 +  1         */
	%       char dim_info;   % char hkey_un0;        /* 39 +  1 */
	%       };                               /* total=40 bytes  */
	%
	% int sizeof_header   Should be 348.
	% char regular        Must be 'r' to indicate that all images and 
	%                     volumes are the same size. 

    v6 = version;
    if str2num(v6(1))<6
       directchar = '*char';
    else
       directchar = 'uchar=>char';
    end

    hk.sizeof_hdr    = fread(fid, 1,'int32')';	% should be 348!
    hk.data_type     = deblank(fread(fid,10,directchar)');
    hk.db_name       = deblank(fread(fid,18,directchar)');
    hk.extents       = fread(fid, 1,'int32')';
    hk.session_error = fread(fid, 1,'int16')';
    hk.regular       = fread(fid, 1,directchar)';
    hk.dim_info      = fread(fid, 1,'uchar')';
    
    return					% header_key


%---------------------------------------------------------------------
function [ dime ] = image_dimension(fid)

	%  Original header structures    
	%  struct image_dimension
	%       {                                /* off + size      */
	%       short int dim[8];                /* 0 + 16          */
        %       /*
        %           dim[0]      Number of dimensions in database; usually 4. 
        %           dim[1]      Image X dimension;  number of *pixels* in an image row. 
        %           dim[2]      Image Y dimension;  number of *pixel rows* in slice. 
        %           dim[3]      Volume Z dimension; number of *slices* in a volume. 
        %           dim[4]      Time points; number of volumes in database
        %       */
	%       float intent_p1;   % char vox_units[4];   /* 16 + 4       */
	%       float intent_p2;   % char cal_units[8];   /* 20 + 4       */
	%       float intent_p3;   % char cal_units[8];   /* 24 + 4       */
	%       short int intent_code;   % short int unused1;   /* 28 + 2 */
	%       short int datatype;              /* 30 + 2          */
	%       short int bitpix;                /* 32 + 2          */
	%       short int slice_start;   % short int dim_un0;   /* 34 + 2 */
	%       float pixdim[8];                 /* 36 + 32         */
	%	/*
	%		pixdim[] specifies the voxel dimensions:
	%		pixdim[1] - voxel width, mm
	%		pixdim[2] - voxel height, mm
	%		pixdim[3] - slice thickness, mm
	%		pixdim[4] - volume timing, in msec
	%					..etc
	%	*/
	%       float vox_offset;                /* 68 + 4          */
	%       float scl_slope;   % float roi_scale;     /* 72 + 4 */
	%       float scl_inter;   % float funused1;      /* 76 + 4 */
	%       short slice_end;   % float funused2;      /* 80 + 2 */
	%       char slice_code;   % float funused2;      /* 82 + 1 */
	%       char xyzt_units;   % float funused2;      /* 83 + 1 */
	%       float cal_max;                   /* 84 + 4          */
	%       float cal_min;                   /* 88 + 4          */
	%       float slice_duration;   % int compressed; /* 92 + 4 */
	%       float toffset;   % int verified;          /* 96 + 4 */
	%       int glmax;                       /* 100 + 4         */
	%       int glmin;                       /* 104 + 4         */
	%       };                               /* total=108 bytes */
	
    dime.dim        = fread(fid,8,'int16')';
    dime.intent_p1  = fread(fid,1,'float32')';
    dime.intent_p2  = fread(fid,1,'float32')';
    dime.intent_p3  = fread(fid,1,'float32')';
    dime.intent_code = fread(fid,1,'int16')';
    dime.datatype   = fread(fid,1,'int16')';
    dime.bitpix     = fread(fid,1,'int16')';
    dime.slice_start = fread(fid,1,'int16')';
    dime.pixdim     = fread(fid,8,'float32')';
    dime.vox_offset = fread(fid,1,'float32')';
    dime.scl_slope  = fread(fid,1,'float32')';
    dime.scl_inter  = fread(fid,1,'float32')';
    dime.slice_end  = fread(fid,1,'int16')';
    dime.slice_code = fread(fid,1,'uchar')';
    dime.xyzt_units = fread(fid,1,'uchar')';
    dime.cal_max    = fread(fid,1,'float32')';
    dime.cal_min    = fread(fid,1,'float32')';
    dime.slice_duration = fread(fid,1,'float32')';
    dime.toffset    = fread(fid,1,'float32')';
    dime.glmax      = fread(fid,1,'int32')';
    dime.glmin      = fread(fid,1,'int32')';
        
    return					% image_dimension


%---------------------------------------------------------------------
function [ hist ] = data_history(fid)
        
	%  Original header structures
	%  struct data_history       
	%       {                                /* off + size      */
	%       char descrip[80];                /* 0 + 80          */
	%       char aux_file[24];               /* 80 + 24         */
	%       short int qform_code;            /* 104 + 2         */
	%       short int sform_code;            /* 106 + 2         */
	%       float quatern_b;                 /* 108 + 4         */
	%       float quatern_c;                 /* 112 + 4         */
	%       float quatern_d;                 /* 116 + 4         */
	%       float qoffset_x;                 /* 120 + 4         */
	%       float qoffset_y;                 /* 124 + 4         */
	%       float qoffset_z;                 /* 128 + 4         */
	%       float srow_x[4];                 /* 132 + 16        */
	%       float srow_y[4];                 /* 148 + 16        */
	%       float srow_z[4];                 /* 164 + 16        */
	%       char intent_name[16];            /* 180 + 16        */
	%       char magic[4];   % int smin;     /* 196 + 4         */
	%       };                               /* total=200 bytes */

    v6 = version;
    if str2num(v6(1))<6
       directchar = '*char';
    else
       directchar = 'uchar=>char';
    end

    hist.descrip     = deblank(fread(fid,80,directchar)');
    hist.aux_file    = deblank(fread(fid,24,directchar)');
    hist.qform_code  = fread(fid,1,'int16')';
    hist.sform_code  = fread(fid,1,'int16')';
    hist.quatern_b   = fread(fid,1,'float32')';
    hist.quatern_c   = fread(fid,1,'float32')';
    hist.quatern_d   = fread(fid,1,'float32')';
    hist.qoffset_x   = fread(fid,1,'float32')';
    hist.qoffset_y   = fread(fid,1,'float32')';
    hist.qoffset_z   = fread(fid,1,'float32')';
    hist.srow_x      = fread(fid,4,'float32')';
    hist.srow_y      = fread(fid,4,'float32')';
    hist.srow_z      = fread(fid,4,'float32')';
    hist.intent_name = deblank(fread(fid,16,directchar)');
    hist.magic       = deblank(fread(fid,4,directchar)');

    fseek(fid,253,'bof');
    hist.originator  = fread(fid, 5,'int16')';
    
    return					% data_history

%  internal function

%  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)

function [img,hdr] = load_nii_img(hdr,filetype,fileprefix,machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB)

   if ~exist('hdr','var') | ~exist('filetype','var') | ~exist('fileprefix','var') | ~exist('machine','var')
      error('Usage: [img,hdr] = load_nii_img(hdr,filetype,fileprefix,machine,[img_idx],[dim5_idx],[dim6_idx],[dim7_idx],[old_RGB]);');
   end

   if ~exist('img_idx','var') | isempty(img_idx) | hdr.dime.dim(5)<1
      img_idx = [];
   end

   if ~exist('dim5_idx','var') | isempty(dim5_idx) | hdr.dime.dim(6)<1
      dim5_idx = [];
   end

   if ~exist('dim6_idx','var') | isempty(dim6_idx) | hdr.dime.dim(7)<1
      dim6_idx = [];
   end

   if ~exist('dim7_idx','var') | isempty(dim7_idx) | hdr.dime.dim(8)<1
      dim7_idx = [];
   end

   if ~exist('old_RGB','var') | isempty(old_RGB)
      old_RGB = 0;
   end

   %  check img_idx
   %
   if ~isempty(img_idx) & ~isnumeric(img_idx)
      error('"img_idx" should be a numerical array.');
   end

   if length(unique(img_idx)) ~= length(img_idx)
      error('Duplicate image index in "img_idx"');
   end

   if ~isempty(img_idx) & (min(img_idx) < 1 | max(img_idx) > hdr.dime.dim(5))
      max_range = hdr.dime.dim(5);

      if max_range == 1
         error(['"img_idx" should be 1.']);
      else
         range = ['1 ' num2str(max_range)];
         error(['"img_idx" should be an integer within the range of [' range '].']);
      end
   end

   %  check dim5_idx
   %
   if ~isempty(dim5_idx) & ~isnumeric(dim5_idx)
      error('"dim5_idx" should be a numerical array.');
   end

   if length(unique(dim5_idx)) ~= length(dim5_idx)
      error('Duplicate index in "dim5_idx"');
   end

   if ~isempty(dim5_idx) & (min(dim5_idx) < 1 | max(dim5_idx) > hdr.dime.dim(6))
      max_range = hdr.dime.dim(6);

      if max_range == 1
         error(['"dim5_idx" should be 1.']);
      else
         range = ['1 ' num2str(max_range)];
         error(['"dim5_idx" should be an integer within the range of [' range '].']);
      end
   end

   %  check dim6_idx
   %
   if ~isempty(dim6_idx) & ~isnumeric(dim6_idx)
      error('"dim6_idx" should be a numerical array.');
   end

   if length(unique(dim6_idx)) ~= length(dim6_idx)
      error('Duplicate index in "dim6_idx"');
   end

   if ~isempty(dim6_idx) & (min(dim6_idx) < 1 | max(dim6_idx) > hdr.dime.dim(7))
      max_range = hdr.dime.dim(7);

      if max_range == 1
         error(['"dim6_idx" should be 1.']);
      else
         range = ['1 ' num2str(max_range)];
         error(['"dim6_idx" should be an integer within the range of [' range '].']);
      end
   end

   %  check dim7_idx
   %
   if ~isempty(dim7_idx) & ~isnumeric(dim7_idx)
      error('"dim7_idx" should be a numerical array.');
   end

   if length(unique(dim7_idx)) ~= length(dim7_idx)
      error('Duplicate index in "dim7_idx"');
   end

   if ~isempty(dim7_idx) & (min(dim7_idx) < 1 | max(dim7_idx) > hdr.dime.dim(8))
      max_range = hdr.dime.dim(8);

      if max_range == 1
         error(['"dim7_idx" should be 1.']);
      else
         range = ['1 ' num2str(max_range)];
         error(['"dim7_idx" should be an integer within the range of [' range '].']);
      end
   end

   [img,hdr] = read_image(hdr,filetype,fileprefix,machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB);

   return					% load_nii_img


%---------------------------------------------------------------------
function [img,hdr] = read_image(hdr,filetype,fileprefix,machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB)

   switch filetype
   case {0, 1}
      fn = [fileprefix '.img'];
   case 2
      fn = [fileprefix '.nii'];
   end

   fid = fopen(fn,'r',machine);

   if fid < 0,
      msg = sprintf('Cannot open file %s.',fn);
      error(msg);
   end

   %  Set bitpix according to datatype
   %
   %  /*Acceptable values for datatype are*/ 
   %
   %     0 None                     (Unknown bit per voxel) % DT_NONE, DT_UNKNOWN 
   %     1 Binary                         (ubit1, bitpix=1) % DT_BINARY 
   %     2 Unsigned char         (uchar or uint8, bitpix=8) % DT_UINT8, NIFTI_TYPE_UINT8 
   %     4 Signed short                  (int16, bitpix=16) % DT_INT16, NIFTI_TYPE_INT16 
   %     8 Signed integer                (int32, bitpix=32) % DT_INT32, NIFTI_TYPE_INT32 
   %    16 Floating point    (single or float32, bitpix=32) % DT_FLOAT32, NIFTI_TYPE_FLOAT32 
   %    32 Complex, 2 float32      (Use float32, bitpix=64) % DT_COMPLEX64, NIFTI_TYPE_COMPLEX64
   %    64 Double precision  (double or float64, bitpix=64) % DT_FLOAT64, NIFTI_TYPE_FLOAT64 
   %   128 uint8 RGB                 (Use uint8, bitpix=24) % DT_RGB24, NIFTI_TYPE_RGB24 
   %   256 Signed char            (schar or int8, bitpix=8) % DT_INT8, NIFTI_TYPE_INT8 
   %   511 Single RGB              (Use float32, bitpix=96) % DT_RGB96, NIFTI_TYPE_RGB96
   %   512 Unsigned short               (uint16, bitpix=16) % DT_UNINT16, NIFTI_TYPE_UNINT16 
   %   768 Unsigned integer             (uint32, bitpix=32) % DT_UNINT32, NIFTI_TYPE_UNINT32 
   %  1024 Signed long long              (int64, bitpix=64) % DT_INT64, NIFTI_TYPE_INT64
   %  1280 Unsigned long long           (uint64, bitpix=64) % DT_UINT64, NIFTI_TYPE_UINT64 
   %  1536 Long double, float128  (Unsupported, bitpix=128) % DT_FLOAT128, NIFTI_TYPE_FLOAT128 
   %  1792 Complex128, 2 float64  (Use float64, bitpix=128) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 
   %  2048 Complex256, 2 float128 (Unsupported, bitpix=256) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128 
   %
   switch hdr.dime.datatype
   case   1,
      hdr.dime.bitpix = 1;  precision = 'ubit1';
   case   2,
      hdr.dime.bitpix = 8;  precision = 'uint8';
   case   4,
      hdr.dime.bitpix = 16; precision = 'int16';
   case   8,
      hdr.dime.bitpix = 32; precision = 'int32';
   case  16,
      hdr.dime.bitpix = 32; precision = 'float32';
   case  32,
      hdr.dime.bitpix = 64; precision = 'float32';
   case  64,
      hdr.dime.bitpix = 64; precision = 'float64';
   case 128,
      hdr.dime.bitpix = 24; precision = 'uint8';
   case 256 
      hdr.dime.bitpix = 8;  precision = 'int8';
   case 511 
      hdr.dime.bitpix = 96; precision = 'float32';
   case 512 
      hdr.dime.bitpix = 16; precision = 'uint16';
   case 768 
      hdr.dime.bitpix = 32; precision = 'uint32';
   case 1024
      hdr.dime.bitpix = 64; precision = 'int64';
   case 1280
      hdr.dime.bitpix = 64; precision = 'uint64';
   case 1792,
      hdr.dime.bitpix = 128; precision = 'float64';
   otherwise
      fclose(fid);
      error('This datatype is not supported'); 
   end

   hdr.dime.dim(find(hdr.dime.dim < 1)) = 1;

   %  move pointer to the start of image block
   %
   switch filetype
   case {0, 1}
      fseek(fid, 0, 'bof');
   case 2
      fseek(fid, hdr.dime.vox_offset, 'bof');
   end

   %  Load whole image block for old Analyze format or binary image;
   %  otherwise, load images that are specified in img_idx, dim5_idx,
   %  dim6_idx, and dim7_idx
   %
   %  For binary image, we have to read all because pos can not be
   %  seeked in bit and can not be calculated the way below.
   %
   if hdr.dime.datatype == 1 | isequal(hdr.dime.dim(5:8),ones(1,4)) | ...
	(isempty(img_idx) & isempty(dim5_idx) & isempty(dim6_idx) & isempty(dim7_idx))

      %  For each frame, precision of value will be read 
      %  in img_siz times, where img_siz is only the 
      %  dimension size of an image, not the byte storage
      %  size of an image.
      %
      img_siz = prod(hdr.dime.dim(2:8));

      %  For complex float32 or complex float64, voxel values
      %  include [real, imag]
      %
      if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
         img_siz = img_siz * 2;
      end
	 
      %MPH: For RGB24, voxel values include 3 separate color planes
      %
      if hdr.dime.datatype == 128 | hdr.dime.datatype == 511
	 img_siz = img_siz * 3;
      end

      img = fread(fid, img_siz, sprintf('*%s',precision));

      d1 = hdr.dime.dim(2);
      d2 = hdr.dime.dim(3);
      d3 = hdr.dime.dim(4);
      d4 = hdr.dime.dim(5);
      d5 = hdr.dime.dim(6);
      d6 = hdr.dime.dim(7);
      d7 = hdr.dime.dim(8);

      if isempty(img_idx)
         img_idx = 1:d4;
      end

      if isempty(dim5_idx)
         dim5_idx = 1:d5;
      end

      if isempty(dim6_idx)
         dim6_idx = 1:d6;
      end

      if isempty(dim7_idx)
         dim7_idx = 1:d7;
      end
   else

      d1 = hdr.dime.dim(2);
      d2 = hdr.dime.dim(3);
      d3 = hdr.dime.dim(4);
      d4 = hdr.dime.dim(5);
      d5 = hdr.dime.dim(6);
      d6 = hdr.dime.dim(7);
      d7 = hdr.dime.dim(8);

      if isempty(img_idx)
         img_idx = 1:d4;
      end

      if isempty(dim5_idx)
         dim5_idx = 1:d5;
      end

      if isempty(dim6_idx)
         dim6_idx = 1:d6;
      end

      if isempty(dim7_idx)
         dim7_idx = 1:d7;
      end

      %  compute size of one image
      %
      img_siz = prod(hdr.dime.dim(2:4));

      %  For complex float32 or complex float64, voxel values
      %  include [real, imag]
      %
      if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
         img_siz = img_siz * 2;
      end

      %MPH: For RGB24, voxel values include 3 separate color planes
      %
      if hdr.dime.datatype == 128 | hdr.dime.datatype == 511
         img_siz = img_siz * 3;
      end

      % preallocate img
      img = zeros(img_siz, length(img_idx)*length(dim5_idx)*length(dim6_idx)*length(dim7_idx) );
      currentIndex = 1;

      for i7=1:length(dim7_idx)
         for i6=1:length(dim6_idx)
            for i5=1:length(dim5_idx)
               for t=1:length(img_idx)

                  %  Position is seeked in bytes. To convert dimension size
                  %  to byte storage size, hdr.dime.bitpix/8 will be
                  %  applied.
                  %
                  pos = sub2ind([d1 d2 d3 d4 d5 d6 d7], 1, 1, 1, ...
			img_idx(t), dim5_idx(i5),dim6_idx(i6),dim7_idx(i7)) -1;
                  pos = pos * hdr.dime.bitpix/8;

                  if filetype == 2
                     fseek(fid, pos + hdr.dime.vox_offset, 'bof');
                  else
                     fseek(fid, pos, 'bof');
                  end

                  %  For each frame, fread will read precision of value
                  %  in img_siz times
                  %
                  img(:,currentIndex) = fread(fid, img_siz, sprintf('*%s',precision));
                  currentIndex = currentIndex +1;

               end
            end
         end
      end
   end

   %  For complex float32 or complex float64, voxel values
   %  include [real, imag]
   %
   if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
      img = reshape(img, [2, length(img)/2]);
      img = complex(img(1,:)', img(2,:)');
   end

   fclose(fid);

   %  Update the global min and max values 
   %
   hdr.dime.glmax = double(max(img(:)));
   hdr.dime.glmin = double(min(img(:)));

   %  old_RGB treat RGB slice by slice, now it is treated voxel by voxel
   %
   if old_RGB & hdr.dime.datatype == 128 & hdr.dime.bitpix == 24
      % remove squeeze
      img = (reshape(img, [hdr.dime.dim(2:3) 3 hdr.dime.dim(4) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
      img = permute(img, [1 2 4 3 5 6 7 8]);
   elseif hdr.dime.datatype == 128 & hdr.dime.bitpix == 24
      % remove squeeze
      img = (reshape(img, [3 hdr.dime.dim(2:4) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
      img = permute(img, [2 3 4 1 5 6 7 8]);
   elseif hdr.dime.datatype == 511 & hdr.dime.bitpix == 96
      img = double(img(:));
      img = single((img - min(img))/(max(img) - min(img)));
      % remove squeeze
      img = (reshape(img, [3 hdr.dime.dim(2:4) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
      img = permute(img, [2 3 4 1 5 6 7 8]);
   else
      % remove squeeze
      img = (reshape(img, [hdr.dime.dim(2:4) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
   end

   if ~isempty(img_idx)
      hdr.dime.dim(5) = length(img_idx);
   end

   if ~isempty(dim5_idx)
      hdr.dime.dim(6) = length(dim5_idx);
   end

   if ~isempty(dim6_idx)
      hdr.dime.dim(7) = length(dim6_idx);
   end

   if ~isempty(dim7_idx)
      hdr.dime.dim(8) = length(dim7_idx);
   end

   return						% read_image

%  internal function

%  'xform_nii.m' is an internal function called by "load_nii.m", so
%  you do not need run this program by yourself. It does simplified
%  NIfTI sform/qform affine transform, and supports some of the 
%  affine transforms, including translation, reflection, and 
%  orthogonal rotation (N*90 degree).
%
%  For other affine transforms, e.g. any degree rotation, shearing
%  etc. you will have to use the included 'reslice_nii.m' program
%  to reslice the image volume. 'reslice_nii.m' is not called by
%  any other program, and you have to run 'reslice_nii.m' explicitly
%  for those NIfTI files that you want to reslice them.
%
%  Since 'xform_nii.m' does not involve any interpolation or any
%  slice change, the original image volume is supposed to be
%  untouched, although it is translated, reflected, or even 
%  orthogonally rotated, based on the affine matrix in the
%  NIfTI header.
%
%  However, the affine matrix in the header of a lot NIfTI files
%  contain slightly non-orthogonal rotation. Therefore, optional
%  input parameter 'tolerance' is used to allow some distortion
%  in the loaded image for any non-orthogonal rotation or shearing
%  of NIfTI affine matrix. If you set 'tolerance' to 0, it means
%  that you do not allow any distortion. If you set 'tolerance' to
%  1, it means that you do not care any distortion. The image will
%  fail to be loaded if it can not be tolerated. The tolerance will
%  be set to 0.1 (10%), if it is default or empty.
%
%  Because 'reslice_nii.m' has to perform 3D interpolation, it can
%  be slow depending on image size and affine matrix in the header.
%  
%  After you perform the affine transform, the 'nii' structure
%  generated from 'xform_nii.m' or new NIfTI file created from
%  'reslice_nii.m' will be in RAS orientation, i.e. X axis from
%  Left to Right, Y axis from Posterior to Anterior, and Z axis
%  from Inferior to Superior.
%
%  NOTE: This function should be called immediately after load_nii.
%  
%  Usage: [ nii ] = xform_nii(nii, [tolerance], [preferredForm])
%  
%  nii	- NIFTI structure (returned from load_nii)
%
%  tolerance (optional) - distortion allowed for non-orthogonal rotation
%	or shearing in NIfTI affine matrix. It will be set to 0.1 (10%),
%	if it is default or empty.
%
%  preferredForm (optional)  -  selects which transformation from voxels
%	to RAS coordinates; values are s,q,S,Q.  Lower case s,q indicate
%	"prefer sform or qform, but use others if preferred not present". 
%	Upper case indicate the program is forced to use the specificied
%	tranform or fail loading.  'preferredForm' will be 's', if it is
%	default or empty.	- Jeff Gunter
%  
%  NIFTI data format can be found on: http://nifti.nimh.nih.gov
%
%  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)
%
function nii = xform_nii(nii, tolerance, preferredForm)

   %  save a copy of the header as it was loaded.  This is the
   %  header before any sform, qform manipulation is done.
   %
   nii.original.hdr = nii.hdr;

   if ~exist('tolerance','var') | isempty(tolerance)
      tolerance = 0.1;
   elseif(tolerance<=0)
      tolerance = eps;
   end

   if ~exist('preferredForm','var') | isempty(preferredForm)
      preferredForm= 's';				% Jeff
   end

   %  if scl_slope field is nonzero, then each voxel value in the
   %  dataset should be scaled as: y = scl_slope * x + scl_inter
   %  I bring it here because hdr will be modified by change_hdr.
   %
   if nii.hdr.dime.scl_slope ~= 0 & ...
	ismember(nii.hdr.dime.datatype, [2,4,8,16,64,256,512,768]) & ...
	(nii.hdr.dime.scl_slope ~= 1 | nii.hdr.dime.scl_inter ~= 0)

      nii.img = ...
	nii.hdr.dime.scl_slope * double(nii.img) + nii.hdr.dime.scl_inter;

      if nii.hdr.dime.datatype == 64

         nii.hdr.dime.datatype = 64;
         nii.hdr.dime.bitpix = 64;
      else
         nii.img = single(nii.img);

         nii.hdr.dime.datatype = 16;
         nii.hdr.dime.bitpix = 32;
      end

      nii.hdr.dime.glmax = max(double(nii.img(:)));
      nii.hdr.dime.glmin = min(double(nii.img(:)));

      %  set scale to non-use, because it is applied in xform_nii
      %
      nii.hdr.dime.scl_slope = 0;

   end

   %  However, the scaling is to be ignored if datatype is DT_RGB24.

   %  If datatype is a complex type, then the scaling is to be applied
   %  to both the real and imaginary parts.
   %
   if nii.hdr.dime.scl_slope ~= 0 & ...
	ismember(nii.hdr.dime.datatype, [32,1792])

      nii.img = ...
	nii.hdr.dime.scl_slope * double(nii.img) + nii.hdr.dime.scl_inter;

      if nii.hdr.dime.datatype == 32
         nii.img = single(nii.img);
      end

      nii.hdr.dime.glmax = max(double(nii.img(:)));
      nii.hdr.dime.glmin = min(double(nii.img(:)));

      %  set scale to non-use, because it is applied in xform_nii
      %
      nii.hdr.dime.scl_slope = 0;

   end

   %  There is no need for this program to transform Analyze data
   %
   if nii.filetype == 0 & exist([nii.fileprefix '.mat'],'file')
      load([nii.fileprefix '.mat']);	% old SPM affine matrix
      R=M(1:3,1:3);
      T=M(1:3,4);
      T=R*ones(3,1)+T;
      M(1:3,4)=T;
      nii.hdr.hist.qform_code=0;
      nii.hdr.hist.sform_code=1;
      nii.hdr.hist.srow_x=M(1,:);
      nii.hdr.hist.srow_y=M(2,:);
      nii.hdr.hist.srow_z=M(3,:);
   elseif nii.filetype == 0
      nii.hdr.hist.rot_orient = [];
      nii.hdr.hist.flip_orient = [];
      return;				% no sform/qform for Analyze format
   end

   hdr = nii.hdr;

   [hdr,orient]=change_hdr(hdr,tolerance,preferredForm);

   %  flip and/or rotate image data
   %
   if ~isequal(orient, [1 2 3])

      old_dim = hdr.dime.dim([2:4]);

      %  More than 1 time frame
      %
      if ndims(nii.img) > 3
         pattern = 1:prod(old_dim);
      else
         pattern = [];
      end

      if ~isempty(pattern)
         pattern = reshape(pattern, old_dim);
      end

      %  calculate for rotation after flip
      %
      rot_orient = mod(orient + 2, 3) + 1;

      %  do flip:
      %
      flip_orient = orient - rot_orient;

      for i = 1:3
         if flip_orient(i)
            if ~isempty(pattern)
               pattern = flipdim(pattern, i);
            else
               nii.img = flipdim(nii.img, i);
            end
         end
      end

      %  get index of orient (rotate inversely)
      %
      [tmp rot_orient] = sort(rot_orient);

      new_dim = old_dim;
      new_dim = new_dim(rot_orient);
      hdr.dime.dim([2:4]) = new_dim;

      new_pixdim = hdr.dime.pixdim([2:4]);
      new_pixdim = new_pixdim(rot_orient);
      hdr.dime.pixdim([2:4]) = new_pixdim;

      %  re-calculate originator
      %
      tmp = hdr.hist.originator([1:3]);
      tmp = tmp(rot_orient);
      flip_orient = flip_orient(rot_orient);

      for i = 1:3
         if flip_orient(i) & ~isequal(tmp(i), 0)
            tmp(i) = new_dim(i) - tmp(i) + 1;
         end
      end

      hdr.hist.originator([1:3]) = tmp;
      hdr.hist.rot_orient = rot_orient;
      hdr.hist.flip_orient = flip_orient;

      %  do rotation:
      %
      if ~isempty(pattern)
         pattern = permute(pattern, rot_orient);
         pattern = pattern(:);

         if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792 | ...
		hdr.dime.datatype == 128 | hdr.dime.datatype == 511

            tmp = reshape(nii.img(:,:,:,1), [prod(new_dim) hdr.dime.dim(5:8)]);
            tmp = tmp(pattern, :);
            nii.img(:,:,:,1) = reshape(tmp, [new_dim       hdr.dime.dim(5:8)]);

            tmp = reshape(nii.img(:,:,:,2), [prod(new_dim) hdr.dime.dim(5:8)]);
            tmp = tmp(pattern, :);
            nii.img(:,:,:,2) = reshape(tmp, [new_dim       hdr.dime.dim(5:8)]);

            if hdr.dime.datatype == 128 | hdr.dime.datatype == 511
               tmp = reshape(nii.img(:,:,:,3), [prod(new_dim) hdr.dime.dim(5:8)]);
               tmp = tmp(pattern, :);
               nii.img(:,:,:,3) = reshape(tmp, [new_dim       hdr.dime.dim(5:8)]);
            end

         else
            nii.img = reshape(nii.img, [prod(new_dim) hdr.dime.dim(5:8)]);
            nii.img = nii.img(pattern, :);
            nii.img = reshape(nii.img, [new_dim       hdr.dime.dim(5:8)]);
         end
      else
         if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792 | ...
		hdr.dime.datatype == 128 | hdr.dime.datatype == 511

            nii.img(:,:,:,1) = permute(nii.img(:,:,:,1), rot_orient);
            nii.img(:,:,:,2) = permute(nii.img(:,:,:,2), rot_orient);

            if hdr.dime.datatype == 128 | hdr.dime.datatype == 511
               nii.img(:,:,:,3) = permute(nii.img(:,:,:,3), rot_orient);
            end
         else
            nii.img = permute(nii.img, rot_orient);
         end
      end
   else
      hdr.hist.rot_orient = [];
      hdr.hist.flip_orient = [];
   end

   nii.hdr = hdr;

   return;					% xform_nii


%-----------------------------------------------------------------------
function [hdr, orient] = change_hdr(hdr, tolerance, preferredForm)

   orient = [1 2 3];
   affine_transform = 1;

   %  NIFTI can have both sform and qform transform. This program
   %  will check sform_code prior to qform_code by default.
   %
   %  If user specifys "preferredForm", user can then choose the
   %  priority.					- Jeff
   %
   useForm=[];					% Jeff

   if isequal(preferredForm,'S')
       if isequal(hdr.hist.sform_code,0)
           error('User requires sform, sform not set in header');
       else
           useForm='s';
       end
   end						% Jeff

   if isequal(preferredForm,'Q')
       if isequal(hdr.hist.qform_code,0)
           error('User requires qform, qform not set in header');
       else
           useForm='q';
       end
   end						% Jeff

   if isequal(preferredForm,'s')
       if hdr.hist.sform_code > 0
           useForm='s';
       elseif hdr.hist.qform_code > 0
           useForm='q';
       end
   end						% Jeff
   
   if isequal(preferredForm,'q')
       if hdr.hist.qform_code > 0
           useForm='q';
       elseif hdr.hist.sform_code > 0
           useForm='s';
       end
   end						% Jeff

   if isequal(useForm,'s')
      R = [hdr.hist.srow_x(1:3)
           hdr.hist.srow_y(1:3)
           hdr.hist.srow_z(1:3)];

      T = [hdr.hist.srow_x(4)
           hdr.hist.srow_y(4)
           hdr.hist.srow_z(4)];

      if det(R) == 0 | ~isequal(R(find(R)), sum(R)')
         hdr.hist.old_affine = [ [R;[0 0 0]] [T;1] ];
         R_sort = sort(abs(R(:)));
         R( find( abs(R) < tolerance*min(R_sort(end-2:end)) ) ) = 0;
         hdr.hist.new_affine = [ [R;[0 0 0]] [T;1] ];

         if det(R) == 0 | ~isequal(R(find(R)), sum(R)')
            msg = [char(10) char(10) '   Non-orthogonal rotation or shearing '];
            msg = [msg 'found inside the affine matrix' char(10)];
            msg = [msg '   in this NIfTI file. You have 3 options:' char(10) char(10)];
            msg = [msg '   1. Using included ''reslice_nii.m'' program to reslice the NIfTI' char(10)];
            msg = [msg '      file. I strongly recommand this, because it will not cause' char(10)];
            msg = [msg '      negative effect, as long as you remember not to do slice' char(10)];
            msg = [msg '      time correction after using ''reslice_nii.m''.' char(10) char(10)];
            msg = [msg '   2. Using included ''load_untouch_nii.m'' program to load image' char(10)];
            msg = [msg '      without applying any affine geometric transformation or' char(10)];
            msg = [msg '      voxel intensity scaling. This is only for people who want' char(10)];
            msg = [msg '      to do some image processing regardless of image orientation' char(10)];
            msg = [msg '      and to save data back with the same NIfTI header.' char(10) char(10)];
            msg = [msg '   3. Increasing the tolerance to allow more distortion in loaded' char(10)];
            msg = [msg '      image, but I don''t suggest this.' char(10) char(10)];
            msg = [msg '   To get help, please type:' char(10) char(10) '   help reslice_nii.m' char(10)];
            msg = [msg '   help load_untouch_nii.m' char(10) '   help load_nii.m'];
            error(msg);
         end
      end

   elseif isequal(useForm,'q')
      b = hdr.hist.quatern_b;
      c = hdr.hist.quatern_c;
      d = hdr.hist.quatern_d;

      if 1.0-(b*b+c*c+d*d) < 0
         if abs(1.0-(b*b+c*c+d*d)) < 1e-5
            a = 0;
         else
            error('Incorrect quaternion values in this NIFTI data.');
         end
      else
         a = sqrt(1.0-(b*b+c*c+d*d));
      end

      qfac = hdr.dime.pixdim(1);
      if qfac==0, qfac = 1; end
      i = hdr.dime.pixdim(2);
      j = hdr.dime.pixdim(3);
      k = qfac * hdr.dime.pixdim(4);

      R = [a*a+b*b-c*c-d*d     2*b*c-2*a*d        2*b*d+2*a*c
           2*b*c+2*a*d         a*a+c*c-b*b-d*d    2*c*d-2*a*b
           2*b*d-2*a*c         2*c*d+2*a*b        a*a+d*d-c*c-b*b];

      T = [hdr.hist.qoffset_x
           hdr.hist.qoffset_y
           hdr.hist.qoffset_z];

      %  qforms are expected to generate rotation matrices R which are
      %  det(R) = 1; we'll make sure that happens.
      %  
      %  now we make the same checks as were done above for sform data
      %  BUT we do it on a transform that is in terms of voxels not mm;
      %  after we figure out the angles and squash them to closest 
      %  rectilinear direction. After that, the voxel sizes are then
      %  added.
      %
      %  This part is modified by Jeff Gunter.
      %
      if det(R) == 0 | ~isequal(R(find(R)), sum(R)')

         %  det(R) == 0 is not a common trigger for this ---
         %  R(find(R)) is a list of non-zero elements in R; if that
         %  is straight (not oblique) then it should be the same as 
         %  columnwise summation. Could just as well have checked the
         %  lengths of R(find(R)) and sum(R)' (which should be 3)
         %
         hdr.hist.old_affine = [ [R * diag([i j k]);[0 0 0]] [T;1] ];
         R_sort = sort(abs(R(:)));
         R( find( abs(R) < tolerance*min(R_sort(end-2:end)) ) ) = 0;
         R = R * diag([i j k]);
         hdr.hist.new_affine = [ [R;[0 0 0]] [T;1] ];

         if det(R) == 0 | ~isequal(R(find(R)), sum(R)')
            msg = [char(10) char(10) '   Non-orthogonal rotation or shearing '];
            msg = [msg 'found inside the affine matrix' char(10)];
            msg = [msg '   in this NIfTI file. You have 3 options:' char(10) char(10)];
            msg = [msg '   1. Using included ''reslice_nii.m'' program to reslice the NIfTI' char(10)];
            msg = [msg '      file. I strongly recommand this, because it will not cause' char(10)];
            msg = [msg '      negative effect, as long as you remember not to do slice' char(10)];
            msg = [msg '      time correction after using ''reslice_nii.m''.' char(10) char(10)];
            msg = [msg '   2. Using included ''load_untouch_nii.m'' program to load image' char(10)];
            msg = [msg '      without applying any affine geometric transformation or' char(10)];
            msg = [msg '      voxel intensity scaling. This is only for people who want' char(10)];
            msg = [msg '      to do some image processing regardless of image orientation' char(10)];
            msg = [msg '      and to save data back with the same NIfTI header.' char(10) char(10)];
            msg = [msg '   3. Increasing the tolerance to allow more distortion in loaded' char(10)];
            msg = [msg '      image, but I don''t suggest this.' char(10) char(10)];
            msg = [msg '   To get help, please type:' char(10) char(10) '   help reslice_nii.m' char(10)];
            msg = [msg '   help load_untouch_nii.m' char(10) '   help load_nii.m'];
            error(msg);
         end

      else
         R = R * diag([i j k]);
      end					% 1st det(R)

   else
      affine_transform = 0;	% no sform or qform transform
   end

   if affine_transform == 1
      voxel_size = abs(sum(R,1));
      inv_R = inv(R);
      originator = inv_R*(-T)+1;
      orient = get_orient(inv_R);

      %  modify pixdim and originator
      %
      hdr.dime.pixdim(2:4) = voxel_size;
      hdr.hist.originator(1:3) = originator;

      %  set sform or qform to non-use, because they have been
      %  applied in xform_nii
      %
      hdr.hist.qform_code = 0;
      hdr.hist.sform_code = 0;
   end

   %  apply space_unit to pixdim if not 1 (mm)
   %
   space_unit = get_units(hdr);

   if space_unit ~= 1
      hdr.dime.pixdim(2:4) = hdr.dime.pixdim(2:4) * space_unit;

      %  set space_unit of xyzt_units to millimeter, because
      %  voxel_size has been re-scaled
      %
      hdr.dime.xyzt_units = char(bitset(hdr.dime.xyzt_units,1,0));
      hdr.dime.xyzt_units = char(bitset(hdr.dime.xyzt_units,2,1));
      hdr.dime.xyzt_units = char(bitset(hdr.dime.xyzt_units,3,0));
   end

   hdr.dime.pixdim = abs(hdr.dime.pixdim);

   return;					% change_hdr


%-----------------------------------------------------------------------
function orient = get_orient(R)

   orient = [];

   for i = 1:3
      switch find(R(i,:)) * sign(sum(R(i,:)))
      case 1
         orient = [orient 1];		% Left to Right
      case 2
         orient = [orient 2];		% Posterior to Anterior
      case 3
         orient = [orient 3];		% Inferior to Superior
      case -1
         orient = [orient 4];		% Right to Left
      case -2
         orient = [orient 5];		% Anterior to Posterior
      case -3
         orient = [orient 6];		% Superior to Inferior
      end
   end

   return;					% get_orient


%-----------------------------------------------------------------------
function [space_unit, time_unit] = get_units(hdr)

   switch bitand(hdr.dime.xyzt_units, 7)	% mask with 0x07
   case 1
      space_unit = 1e+3;		% meter, m
   case 3
      space_unit = 1e-3;		% micrometer, um
   otherwise
      space_unit = 1;			% millimeter, mm
   end

   switch bitand(hdr.dime.xyzt_units, 56)	% mask with 0x38
   case 16
      time_unit = 1e-3;			% millisecond, ms
   case 24
      time_unit = 1e-6;			% microsecond, us
   otherwise
      time_unit = 1;			% second, s
   end

   return;					% get_units

