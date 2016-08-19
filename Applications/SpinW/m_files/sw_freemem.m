function mem = sw_freemem()
% gives the amount of free RAM in bytes
%
% mem = SW_FREEMEM()
%
% If it cannot determine the size of the free memory, it returns zero.
%

% $Name: SpinW$ ($Version: 2.1$)
% $Author: S. Toth$ ($Contact: sandor.toth@psi.ch$)
% $Revision: 238 $ ($Date: 07-Feb-2015 $)
% $License: GNU GENERAL PUBLIC LICENSE$

mem = 0;

try %#ok<TRYNC>
    if ispc
        memStr = memory;
        mem = memStr.MemAvailableAllArrays;
    elseif ismac
        [~,memStr] = unix('vm_stat | grep free');
        mem = sscanf(memStr(14:end),'%f')*4096;
    elseif isunix
        [~, memStr] = unix('free -b | grep ''-''');
        [~, mem_free] = strtok(memStr(20:end));
        mem = str2double(mem_free);
    end
end

end