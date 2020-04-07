function s = read_hspy(filename)
% READ_HSPY Read an HypersPy HDF5 data set.
%   See https://hyperspy.org/

  s = [];

  if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
    hspy.name        = 'HypersPy HDF5 data format';
    hspy.method      = mfilename;
    hspy.extension   ={'hspy'};
    s = hspy;
    return
  end
  
  [s.data, s.axes] = readHyperSpyH5(filename);

% ------------------------------------------------------------------------------
% https://github.com/jat255/readHyperSpyH5
% 2017 Joshua Taillon
%
% Terms of Use
% This software was developed at the National Institute of Standards and Technology
% by employees of the Federal Government in the course of their official duties.
% Pursuant to title 17 section 105 of the United States Code this software is not
% subject to copyright protection and is in the public domain. readHyperSpyH5 is
% an experimental software, and NIST assumes no responsibility whatsoever for its
% use by other parties, and makes no guarantees, expressed or implied, about its
% quality, reliability, or any other characteristic. We would appreciate
% acknowledgement if the software is used.
%
% This software can be redistributed and/or modified freely provided that any
% derivative works bear some notice that they are derived from it, and any
% modified versions bear some notice that they have been modified.

function [data, ax] = readHyperSpyH5(filename)
%readHyperSpyH5 Read an HDF5 file saved by HyperSpy
%------------------------------------------------------------------
% Input:
%   filename: name of HyperSpy HDF5 file to read
%   
% Output:
%   data: the numerical data contained within the file
%   ax.scale: double array containing the scales of each axis
%   ax.units: cell array containing the units of each axis
%   ax.name: cell array containing the names of each axis
%   ax.size: double array indicated length of each axis
%   ax.offset:  double array indicating axis offset value
%   ax.navigate: logical array to show if each axis was a nav axis or not

  persistent h5_present % EF 2020: allow faster execution.

  if isempty(h5_present)
    if exist('h5info')
      h5_present = 1;
    else
      h5_present = 0;
    end
  end

  % open low-level h5 file access
  fid = H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');

  % get name of dataset
  if h5_present
    info = h5info(filename, '/Experiments');
  else
    info = hdf5info(filename);  % EF 2020: compatibility with old Matlab.
    info = info.GroupHierarchy;
    % get Experiments group
    if strcmp(info.Name, '/') && strcmp(info.Groups.Name, '/Experiments')
      info = info.Groups;
    end
  end
  name = info.Groups.Name;

  if isfield(info.Groups.Datasets,'ChunkSize')
    numAxes = size(info.Groups.Datasets.ChunkSize,2);
  else
    numAxes = size(info.Groups.Datasets.Chunksize,2);
  end

  ax = [];

  ax.scale    = zeros(1, numAxes);
  ax.units    = cell(1, numAxes);
  ax.name     = cell(1, numAxes);
  ax.size     = zeros(1, numAxes);
  ax.offset   = zeros(1, numAxes);
  ax.navigate = false(1, numAxes);

  % for a 1D [x, y, z, signal] HyperSpy signal, these values all read in
  % as [z, y, x, signal]
  % EF 2020: we use a loop to handle all attributes. Cleaner implementation.
  for i = 1:numAxes
    for att={'scale','units','name','size','offset','navigate'}
      this = ax.(att{1});
      if ~exist('h5readatt'), h5readatt = @h5attget; end % use HDF5TOOLS
      val = h5readatt(filename, [name '/axis-' num2str(i-1)], att{1});
      if ~isempty(val)
        if iscell(val) && isscalar(val), val=val{1}; end
        if iscell(this), this{i} = val;
        else this(i) = val; end
      end
      ax.(att{1}) = this;
    end
  end

  % shift the arrays so they're in [x, y, z, ..., signal] order
  % EF 2020: we use a loop to handle all attributes. Cleaner implementation.
  for att={'scale','units','name','size','offset','navigate'}
    ax.(att{1}) = circshift(ax.(att{1}) ,[0,-1]);
  end


  % read actual data
  % for a 1D [x, y, z, signal] HyperSpy signal, this data is read 
  % as shape [signal, x, y, z]
  if ~exist('h5read')
    data = myh5read(filename, ['/data'], name);
  else
    data = h5read(filename, [name '/data']);
  end

  % shift dimensions of data to put in [x, y, z, signal] order (like in HyperSpy)
  data = shiftdim(data, 1);

  % close low-level file access
  H5F.close(fid);

% ------------------------------------------------------------------------------
function [val, fileID] = myh5read(filename, name, root)
% MYH5READ a compatibility function for h5read

  % hdf5read can stall in R2010a. We use a low-level read
  % val = hdf5read(filename,[data_info.Datasets(i).Name]);
  if nargin < 3, root = ''; end
  if isempty(root), root='/'; end
  if ischar(filename)
    fileID = H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
  else
    fileID = filename;
  end
  dataID = [];
  try
    dataID = H5D.open(fileID, [ root name ]);
  end
  try
    if isempty(dataID), dataID = H5D.open(fileID, name); end
  end
  if ~isempty(dataID)
    val    = H5D.read(dataID, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');
    H5D.close(dataID);
  else
    warning([ mfilename ': ' filename ': ignoring invalid DataSet ' name ]);
    val = [];
  end
  % H5F.close(fileID); done after reading all groups

% ------------------------------------------------------------------------------
% https://fr.mathworks.com/matlabcentral/fileexchange/17172-hdf5tools
% John (2020). HDF5TOOLS MATLAB Central File Exchange. Retrieved April 7, 2020.
%
%  Copyright (c) 2016, The MathWorks, Inc.
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are
%  met:
%
%      * Redistributions of source code must retain the above copyright
%        notice, this list of conditions and the following disclaimer.
%      * Redistributions in binary form must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in
%        the documentation and/or other materials provided with the distribution.
%      * In all cases, the software is, and all modifications and derivatives
%        of the software shall be, licensed to you solely for use in conjunction
%        with MathWorks products and service offerings.
%
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%  POSSIBILITY OF SUCH DAMAGE.

function attval = h5attget(hfile,varname,attname)
%H5ATTGET  Retrieve attribute from HDF5 file.
% 
%   ATTVAL = H5ATTGET(HFILE,LOCATION,ATTNAME) retrieves an attribute named
%   ATTNAME from the HDF5 file HFILE.  LOCATION can be either a group
%   or a traditional dataset.  Use '/' for the root group.
%
%   Scalar string attributes are reported as a row string, rather than a
%   column string.  Other string attributes are returned as-is.
% 
%   1D attributes are returned as a single row, rather than as a column.
% 
%   Example:  
%       d = h5attget('example.h5','/g1/g1.1/dset1.1.1','att1');
%
%   See also h5attput.

%   Copyright 2009-2010 The MathWorks, Inc.

parent_is_dataset = true;

%
% Just use the defaults for now?
flags = 'H5F_ACC_RDONLY';
plist_id = 'H5P_DEFAULT';

file_id = H5F.open ( hfile, flags, plist_id );

if strcmp(varname,'/')
    parent_is_dataset = false;
    parent_id = H5G.open ( file_id, varname );
else
    try
        parent_id=H5D.open(file_id,varname);
    catch %#ok<CTCH>
        parent_is_dataset = false;
        parent_id = H5G.open ( file_id, varname );
    end
end

attr_id = H5A.open_name ( parent_id, attname );

attval = H5A.read ( attr_id, 'H5ML_DEFAULT' );

% Read the datatype information and use that to possibly post-process
% the attribute data.
attr_datatype = H5A.get_type(attr_id);
attr_class = H5T.get_class(attr_datatype);

attr_space_id = H5A.get_space(attr_id);
ndims = H5S.get_simple_extent_dims(attr_space_id);

% Do we alter the size, shape of the attribute?
switch (attr_class)

	case H5ML.get_constant_value('H5T_STRING')
		% Is it a scalar attribute?
		if ( ndims == 0 )
			% Yep, it's a scalar.  Transpose it.
			attval = attval';
		end

	case { H5ML.get_constant_value('H5T_INTEGER'), ...
	       H5ML.get_constant_value('H5T_FLOAT') }
		if ( ndims == 1 )
			% If 1D, make it a row vector.
			attval = reshape(attval,1,numel(attval));
		end
end

H5T.close(attr_datatype);


H5A.close(attr_id);

if parent_is_dataset
    H5D.close(parent_id);    
else
    H5G.close(parent_id);    
end
H5F.close(file_id);
