function [b,sigma] = cumsum(a,dim)
% CUMSUM Cumulative sum of elements.
%   S = CUMSUM(A) computes the cumulative sum of objects elements along
%   columns. S has the same size as A.
%
%   CUMSUM(A,DIM) operates along axis of rank DIM.
%
% Example: a=estruct(peaks); s=cumsum(a); abs(sum(s,0)+1.79e3) < 0.5
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plus, estruct/sum, estruct/prod, estruct/cumprod


if nargin < 2, dim=1; end

[b,sigma] = private_sumtrapzproj(a,dim, 'cumsum');
