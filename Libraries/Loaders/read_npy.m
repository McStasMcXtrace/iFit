function s = read_npy(filename)
% data = read_npy(filename)
%
% Read a Python Numpy array '.npy' file.
%   can handle: 'uint8','uint16','uint32','uint64
%               'int8','int16','int32','int64'
%               'single','double', 'logical'
%
% To generate a NPY file, use python:
% >>> import numpy as np
% >>> x = np.arange(10)
% >>> np.save('outfile.npy', x)
%
% References:
% format specification: https://docs.scipy.org/doc/numpy/neps/npy-format.html
% contrib: https://github.com/kwikteam/npy-matlab, C. Rossant 2016


% Function to read NPY files into matlab. 
% *** Only reads a subset of all possible NPY files, specifically N-D arrays of certain data types. 
% See https://github.com/kwikteam/npy-matlab/blob/master/npy.ipynb for
% more. 
%
s = [];

% in private dir.
[s.shape, s.dataType, s.fortranOrder, s.littleEndian, s.totalHeaderLength, s.npyVersion] = readNPYheader(filename);
if isempty(s.shape), s=[]; return; end

if s.littleEndian
    fid = fopen(filename, 'r', 'l');
else
    fid = fopen(filename, 'r', 'b');
end

try
    
    [~] = fread(fid, s.totalHeaderLength, 'uint8');

    % read the data
    data = fread(fid, prod(s.shape), [s.dataType '=>' s.dataType]);
    
    if length(s.shape)>1 && ~s.fortranOrder
        data = reshape(data, s.shape(end:-1:1));
        data = permute(data, [length(s.shape):-1:1]);
    elseif length(s.shape)>1
        data = reshape(data, s.shape);
    end
    s.data = data;
    fclose(fid);
    
catch me
    s = [];
    fclose(fid);
    rethrow(me);
end


