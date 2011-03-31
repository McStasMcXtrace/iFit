function a=load_fig(a0)
% function a=load_fig(a0)
%
% Returns an iData style dataset from a Matlab figure
%

a = iData(a0.Data.Handle);
close(a0.Data.Handle);
