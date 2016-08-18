% Save a sqw object or array of sqw objects to file
%
%   >> save (w)              % prompt for file
%   >> save (w, file)        % give file
%
% Input:
%   w       sqw object
%   file    [optional] File for output. if none given, then prompted for a file
%   
%   Note that if w is an array of sqw objects then file must be a cell
%   array of filenames of the same size.
%
% Output:
%%   Overloaded methods:
%      sqw/save
%      testsigvar/save
%      sigvar/save
%      IX_dataset_3d/save
%      IX_dataset_2d/save
%      IX_dataset_1d/save
%      IX_axis/save
%      IX_mask/save
%      IX_map/save
%      IX_sample/save
%      IX_moderator/save
%      IX_fermi_chopper/save
%      IX_aperture/save
%      spe/save
%      slice/save
%      phxObject/save
%      parObject/save
%      cut/save
%      sqw/save
%      d4d/save
%      d3d/save
%      d2d/save
%      d1d/save
%      d0d/save
%      COM/save
%%   Reference page in Help browser
%      doc save
%