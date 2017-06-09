function out = opennpy(filename)
%OPENEDF Open an python Numpy NPY Format file, display it
%        and set the 'ans' variable to an iData object with its content
%
% To generate a NPY file, use python:
% >>> import numpy as np
% >>> x = np.arange(10)
% >>> np.save('outfile.npy', x)
%
% (c) E.Farhi, ILL. License: EUPL.


if ~isa(filename,'iData')
  out = iData(filename,'NPY');
else
  out = filename;
end
clear filename;

if numel(out) > 1
  % handle input iData arrays
  for index=1:numel(out)
    out(index) = feval(mfilename, out(index));
  end
end

if ~nargout
  figure; subplot(out);
  
  if ~isdeployed
    assignin('base','ans',out);
    ans = out
  end
end

