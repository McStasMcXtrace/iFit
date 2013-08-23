function a = loadobj(a)
% b = loadobj(s) : de-serialize the object for a much faster load/save.
%
%   @iData/loadobj function to de-serialize object
%     which is re-built from its unit8 representation when loading from e.g. MAT.
%
% input:  s: object or array (iData/uint8)
% output: b: de-serialized object (iData)
% ex:     b=loadobj(a);
%
% Version: $Revision: 1035 $
% See also iData, iData/save, iData/saveas, iData/load

% handle input iData arrays
if numel(a) > 1
  parfor index=1:numel(a)
    a(index) = feval(mfilename, a(index));
  end
  return
end

a.Data = hlp_deserialize(a.Data);

