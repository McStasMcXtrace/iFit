% Takes a section out of an sqw object
%
% Syntax:
%   >> wout = section (win, [ax_1_lo, ax_1_hi], [ax_2_lo, ax_2_hi], ...)
%
% Input:
% ------
%   win                 Input sqw object
%
%   [ax_1_lo, ax_1_hi]  Lower and upper limits for the first axis. Bins are retained whose
%                      centres lie in this range.
%                       To retain the limits of the input structure, type '', [], or the scalar '0'
%
%   [ax_2_lo, ax_2_hi]  Lower and upper limits for the second axis
%
%           :                       :
%  
%       for as many axes as there are plot axes
%
% Output:
% -------
%   wout                Output dataset.
%
%
% Example: to alter the limits of the first and third axes of a 3D sqw object:
%   >> wout = section (win, [1.9,2.1], 0, [-0.55,-0.45])
%                                                           
%%   Overloaded methods:
%      sqw/section
%      sqw/section
%      d4d/section
%      d3d/section
%      d2d/section
%      d1d/section
%      d0d/section
%