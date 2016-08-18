% Read sqw object from a file or array of sqw objects from a set of files
% 
%   >> w=read(sqw,file)
%
% Need to give first argument as an sqw object to enforce a call to this function.
% Can simply create a dummy object with a call to sqw:
%    e.g. >> w = read(sqw,'c:\temp\my_file.sqw')
%
% Input:
% -----
%   sqw         Dummy sqw object to enforce the execution of this method.
%               Can simply create a dummy object with a call to sqw:
%                   e.g. >> w = read(sqw,'c:\temp\my_file.sqw')
%
%   file        File name, or cell array of file names. In this case, reads
%               into an array of sqw objects
%
% Output:
% -------
%   w           sqw object, or array of sqw objects if given cell array of
%               file names
%%   Overloaded methods:
%      VideoReader/read
%      sqw/read
%      testsigvar/read
%      sigvar/read
%      IX_dataset_3d/read
%      IX_dataset_2d/read
%      IX_dataset_1d/read
%      IX_axis/read
%      IX_mask/read
%      IX_map/read
%      IX_sample/read
%      IX_moderator/read
%      IX_fermi_chopper/read
%      IX_aperture/read
%      spe/read
%      slice/read
%      phxObject/read
%      parObject/read
%      cut/read
%      sqw/read
%      d4d/read
%      d3d/read
%      d2d/read
%      d1d/read
%      d0d/read
%      mmreader/read
%%   Reference page in Help browser
%      doc read
%