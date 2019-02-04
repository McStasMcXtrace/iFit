function memStats = memoryInfo()
    % Cross-platform function to get memory usage
    % 
    % @author Sviatoslav Danylenko <dev@udf.su>
    % @license BSD
    %
    % Return values:
    % memStats: memory info @type struct.
    %
    % The following is the fields of memStats:
    % free: the amount of free memory in kB
    % swap: the struct with info for swap (paging) file in kB
    %   swap.usedMatlab: the total amount of by current MATLAB process used swap
    %   swap.used: the total amount of used swap
    %   swap.total: the total amount of used swap
    %   swap.free: the total amount of used swap
    % total: the total amount of memory in kB
    % used: the total amount of used memory in kB
    % cache: the total amount of memory used for cache in kB
    % usedMatlab: the maount of by current MATLAB process used memory in kB
    
    % Commands List for Unix:
    %  Get process memory usage:
    %   ps -p <PID> -o rss --no-headers
    %  Get process memory and swap usage in separate rows:
    %   awk '/VmSwap|VmRSS/{print $2}' /proc/<PID>/status
    %  Get memory usage information:
    %   free -m | sed -n ''2p''
    
    % <https://fr.mathworks.com/matlabcentral/fileexchange/45598-memoryinfo>
%Copyright (c) 2014, Sviatoslav Danylenko
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are
%met:
%
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in
%      the documentation and/or other materials provided with the distribution
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%POSSIBILITY OF SUCH DAMAGE.
    
    memStats = struct('total', false, ...
        'free', false, ...
        'used', false, ...
        'cache', false, ...
        'usedMatlab', false, ...
        'swap', struct('total', false, ...
            'free', false, ...
            'used', false, ...
            'usedMatlab', false)...
        );
    if isunix
        pid = feature('getpid');
        
        [~, usedMatlab] = unix(['awk ''/VmSwap|VmRSS/{print $2}'' /proc/' num2str(pid) '/status']);
        try
            usedMatlab = regexp(usedMatlab, '\n', 'split');
            memStats.usedMatlab = str2double(usedMatlab{1});
            memStats.swap.usedMatlab = str2double(usedMatlab{2});
        catch err %#ok
        end
        
        [~, memUsageStr] = unix('free -k | sed -n ''2p''');
        try
            memUsage = cell2mat(textscan(memUsageStr,'%*s %u %u %u %*u %*u %u','delimiter',' ','collectoutput',true,'multipleDelimsAsOne',true));
            memStats.total = memUsage(1);
            memStats.used = memUsage(2);
            memStats.free = memUsage(3);
            memStats.cache = memUsage(4);
        catch err %#ok
        end
        
        
        [~, swapUsageStr] = unix('free -k | tail -1');
        try
            swapUsage = cell2mat(textscan(swapUsageStr,'%*s %u %u %u','delimiter',' ','collectoutput',true,'multipleDelimsAsOne',true));
            memStats.swap.total = swapUsage(1);
            memStats.swap.used = swapUsage(2);
            memStats.swap.free = swapUsage(3);
        catch err %#ok
        end
    else
        [user, sys] = memory;
        memStats.usedMatlab = bytes2kBytes(user.MemUsedMATLAB);
        memStats.total = bytes2kBytes(sys.PhysicalMemory.Total);
        memStats.free = bytes2kBytes(sys.PhysicalMemory.Available);
        memStats.used = memStats.total - memStats.free;
        freeMemWithSwap = bytes2kBytes(sys.SystemMemory.Available);
        freeSwap = freeMemWithSwap - memStats.free;
        memStats.swap.free = freeSwap; % swap available for MATLAB
%         memStats.swap.free = bytes2kBytes(sys.VirtualAddressSpace.Available);
        memStats.swap.total = bytes2kBytes(sys.VirtualAddressSpace.Total);
        memStats.swap.used = memStats.swap.total - memStats.swap.free;
    end
    
end

function kbytes = bytes2kBytes(bytes)
    % bytes to kB konversion
    kbytes = round(bytes/1024);
end
