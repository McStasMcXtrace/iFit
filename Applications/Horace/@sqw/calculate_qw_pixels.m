% Calculate qh,qk,ql,en for the pixels in an sqw dataset
%
%   >> qw=calculate_qw_pixels(win)
%
%   win     Input sqw object
%
%   qw      Components of momentum (in rlu) and energy for each pixel in the dataset
%           Arrays are packaged as cell array of column vectors for convenience
%           with fitting routines etc.
%               i.e. qw{1}=qh, qw{2}=qk, qw{3}=ql, qw{4}=en
%
%%   Overloaded methods:
%      sqw/calculate_qw_pixels
%      sqw/calculate_qw_pixels
%