function a=load_fig(a0)
% function a=load_fig(a0)
%
% Returns an iData style dataset from a Matlab figure
%
% Version: $Date$ $Version$ $Author$
% See also: iData/load, iLoad, save, iData/saveas
% 

if ~isa(a0,'iData')
  a =iData(iLoad(a0,'fig')); % no post-processing
end

% handle input iData arrays
if numel(a0) > 1
  for index=1:numel(a0)
    a(index) = feval(mfilename, a0(index));
  end
  return
end

a = iData(a0.Data.Handle);
close(a0.Data.Handle);
