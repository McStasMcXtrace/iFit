function b = cumtrapz(a,dim)
% CUMTRAPZ Cumulative trapezoidal numerical integration.
%   S = CUMTRAPZ(A) computes the cumulative integral (primitive) of the data set.
%
%   CUMTRAPZ(A,DIM) operates along axis of rank DIM. The axis is then removed.
%   Default is to use DIM=1. If DIM=0, integration is done on all axes and 
%   the total is returned as a scalar value. 
%   CUMTRAPZ is complementary to CUMSUM and CAMPROJ, and takes into account axes.
%
% Example: a=estruct(peaks); t=cumtrapz(a); abs(sum(t,0)+2.1330e+03) < 0.1
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/plus, estruct/sum, estruct/prod, estruct/cumsum

if nargin < 2, dim=1; end

b = private_sumtrapzproj(a,dim, 'cumtrapz');
