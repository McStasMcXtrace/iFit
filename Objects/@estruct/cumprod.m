function [b,sigma] = cumprod(a,dim)
% CUMPROD Cumulative product of elements.
%   P = CUMPROD(A) computes the cumulative product of objects elements along
%   columns. P has the same size as A.
%
%   CUMPROD(A,DIM) operates along axis of rank DIM.
%
% Example: a=estruct(peaks); p=cumprod(a); abs(sum(p,0)-9.1739e+16) < 5e11
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plus, estruct/sum, estruct/prod, estruct/cumprod

if nargin < 2, dim=1; end

[b,sigma] = private_sumtrapzproj(a,dim, 'cumprod');
