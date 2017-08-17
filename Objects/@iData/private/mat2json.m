function J=mat2json(M,F)

%JSON2MAT converts a Matlab structure into a javascript data object (JSON).
%         M can also be a file name. In the spirit of fast prototyping 
%         this function takes a very loose approach to data types and 
%         dimensionality - neither is explicitly retained.
%
%         The second input argument is optional and when used it indicates
%         the name of the file where J is to be stored.
%
%Example: mat2json(json2mat('{lala:2,lele:4,lili:[1,2,{bubu:5}]}')) 
%
% Jonas Almeida, March 2010

%  Copyright (c) 2010, Jonas Almeida
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are
%  met:
%
%      * Redistributions of source code must retain the above copyright
%        notice, this list of conditions and the following disclaimer.
%      * Redistributions in binary form must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in
%        the documentation and/or other materials provided with the distribution
%
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%  POSSIBILITY OF SUCH DAMAGE.

switch class(M)
    case 'struct'
        J='{';
        f=fieldnames(M);
        for i=1:length(f)
            try
              J=[J,'"',f{i},'":',mat2json(M.(f{i})),','];
            catch ME
              disp(ME.message)
              disp([ mfilename ': Ignoring field "' f{i} '" which is not a base Matlab class, but is ' class(M.(f{i})) ]);
            end
        end
        J(end)='}';
        
    case 'cell'
        J='[';
        for i=1:length(M)
            try
              J=[J,mat2json(M{i}),','];
            catch
              disp([ mfilename ': Ignoring cell index "' num2str(i) '" which is not a base Matlab class, but is ' class(M{i}) ]);
            end
        end
        J(end)=']';
    otherwise
        if isnumeric(M) % notice looseness in not converting single numbers into arrays
            if length(M(:))==1
                J=num2str(M);
            else
                s=size(M);
                if numel(s) == 2 && any(s == 1) % vector
                    J=['[', val2str(M)  ,']']; % and of destroying dimensionality
                    J=regexprep(J,'\s+',',');
                elseif length(s)==2 %2D solution
                    J='[';
                    for i=1:s(1)
                        J=[J,mat2json(M(i,:)),','];
                    end
                    J(end)=']';
                elseif length(s)>2 % for now treat higher dimensions as linear vectors
                    J=['[', val2str(M(:)') ,']']; % and of destroying dimensionality
                    J=regexprep(J,'\s+',',');
                end
            end
        else
          try
            J=mat2json(struct(M));  % OK for objects
          catch
            J=['"',char(M),'"']; % otherwise it is treated as a string
          end
        end
end

if nargin>1 %save JSON result in file
    fid=fopen(F,'w');
    fprintf(fid,'%s',J);
    fclose(fid);
end

% ------------------------------------------------------------------------------
function val = val2str(val)
  s = size(val);
  if numel(s) == 2 && s(1) == 1     % vector: horizontal
    val = sprintf('%g ',val);
  elseif numel(s) == 2 && s(2) == 1 % vector: vertical
    val = sprintf('%g\n',val);
  else
    val = num2str(val);
  end
  
          
