function a = saveobj(a)
% b = saveobj(s) : serialize the object for a much faster load/save.
%
%   @iData/saveobj function to serialize object
%     which is converted to unit8 before saving.
%
% input:  s: object or array (iData)
% output: b: serialized object (uint8)
% ex:     b=saveobj(a);
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

a.Data = hlp_serialize(a.Data);

