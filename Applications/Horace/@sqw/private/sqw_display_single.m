% Display useful information from an sqw object
%
% Syntax:
%
%   >> sqw_display_single (din)
%   >> sqw_display_single (din,npixtot,type)
%
%   din             Structure from sqw object (sqw-type or dnd-type)
%
% Optionally:
%   npixtot         total number of pixels if sqw type
%   type            data type: 'a' or 'b+'
%                  
%   If the optional parameters are given, then only the header information
%   part of data needs to be passed, namely the fields:
%      uoffset,u_to_rlu,ulen,ulabel,iax,iint,pax,p,dax[,urange]
%  (urange is only present if sqw type object)
%
%   If the optional parameters are not given, then the whole data structure
%   needs to be given, and npixtot and type are computed from the structure.
%
%   If an optional parameter is given but is empty, then the missing value for that
%   parameter is computed from the data structure.
%
%