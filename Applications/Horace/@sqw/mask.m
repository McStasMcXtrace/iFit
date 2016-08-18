% Remove the bins indicated by the mask array
%
% Syntax:
%   >> wout = mask (win, mask_array)
%
% Input:
% ------
%   win                 Input sqw object
%
%   mask_array          Array of 1 or 0 (or true or false) that indicate
%                      which points to retain (true to retain, false to ignore)
%                       Numeric or logical array of same number of elements
%                      as the data.
%                       Note: mask will be applied to the stored data array
%                      according as the projection axes, not the display axes.
%                      Thus permuting the display axes does not alter the
%                      effect of masking the data.
%
% Output:
% -------
%   wout                Output dataset.
%%   Overloaded methods:
%      sqw/mask
%      IX_dataset_3d/mask
%      IX_dataset_2d/mask
%      IX_dataset_1d/mask
%      sqw/mask
%      d4d/mask
%      d3d/mask
%      d2d/mask
%      d1d/mask
%      d0d/mask
%