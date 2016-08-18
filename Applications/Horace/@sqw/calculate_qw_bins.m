% Calculate qh,qk,ql,en for the centres of the bins of an n-dimensional sqw dataset
%
%   >> qw=calculate_qw_bins(win)
%
% Input:
% ------
%   win     Input sqw object
%
% Output:
% -------
%   qw      Components of momentum (in rlu) and energy for each bin in the dataset
%           Arrays are packaged as cell array of column vectors for convenience
%           with fitting routines etc.
%               i.e. qw{1}=qh, qw{2}=qk, qw{3}=ql, qw{4}=en
%%   Overloaded methods:
%      sqw/calculate_qw_bins
%      sqw/calculate_qw_bins
%      d4d/calculate_qw_bins
%      d3d/calculate_qw_bins
%      d2d/calculate_qw_bins
%      d1d/calculate_qw_bins
%      d0d/calculate_qw_bins
%