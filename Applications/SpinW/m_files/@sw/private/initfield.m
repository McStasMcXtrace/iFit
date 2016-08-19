function obj = initfield(obj)
% SW_INITFIELD(objS) initializes all subfields of the obj structure to
% their initial values.
%

% $Name: SpinW$ ($Version: 2.1$)
% $Author: S. Toth$ ($Contact: sandor.toth@psi.ch$)
% $Revision: 238 $ ($Date: 07-Feb-2015 $)
% $License: GNU GENERAL PUBLIC LICENSE$

datstruct = datastruct();
mainfield = datstruct.mainfield;
subfield  = datstruct.subfield;
defval    = datstruct.defval;


for ii = 1:length(mainfield)
    for jj = 1:size(subfield,2)
        if ~isempty(subfield{ii,jj})
            if ~isfield(obj,mainfield{ii}) || ~isfield(eval(['obj.' mainfield{ii}]),subfield{ii,jj})
                obj.(mainfield{ii}).(subfield{ii,jj}) = eval(defval{ii,jj});
            end
        end
    end
end

end