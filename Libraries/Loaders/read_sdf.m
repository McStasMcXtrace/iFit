% data = read_sdf(filename) Read a HP/Agilent/Keysight Standard Data Format (SDF)
%
% Description:
% Macro for reading HP/Agilent/Keysight Standard Data Format (SDF)
%
% References:
%   27 May 2018 Justin Dinale
%   <https://fr.mathworks.com/matlabcentral/fileexchange/67513-sdf-importer>
% See also: read_edf, read_adsc, read_edf, read_sif, read_cbf, read_spe, read_fits, read_hbin, read_image

function data = read_sdf(filename)

  data = SDF_import(filename);
  
end

% ------------------------------------------------------------------------------
function [data, file_raw] = SDF_import( filename )
%SDF_IMPORT - Imports HP/Agilent/Keysight Standard Data Format (SDF) files to MATLAB and OCTAVE
%  This code is an implementation of the SDF standard as defined in:
%  'User's Guide - Standard Data Format Utilities, Version B.02.01'
%  by Agilent Technologies
%  Manufacturing Part Number: 5963-1715
%  Printed in USA
%  December 1994
%  © Copyright 1989, 1991-94, 2005 Agilent Technologies, Inc.
%
% Syntax: [data file_raw] = SDF_import( filename )
%
% Inputs:
%   filename - String of the filename of the SDF file including the file extension.
%
% Outputs:
%   data - Extracted SDF headers, and processed signal traces.
%   file_raw - Raw byte array of the SDF file which is being intepreted.
%
% Other m-files required: SDF_template.m
% Subfunctions: capitaliseFirstLetters, dataTypeByteSize, SDF_DataType,
%   SDF_uintX, SDF_int8, SDF_char, SDF_short, SDF_long, SDF_float, SDF_double, SDF_struct,
%   SDF_uintX_vector, SDF_short_vector, SDF_long_vector, SDF_float_vector, SDF_double_vector,
%   SDF_DATA_SizesX, SDF_DATA_SizesY, SDF_extract_labels, SDF_Interpret_Template,
%   SDF_Multi_ChannelX, SDF_Multi_ChannelY, SDF_Multi_TraceX, SDF_Multi_TraceY,
%   SDF_Single_TraceX, SDF_Single_TraceY, SDF_XDATA_Process, SDF_YDATA_Process, SDF_Y_scale.
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 27/04/2018
% © Copyright 2016-2018 Commonwealth of Australia, represented by the Department of Defence
%
% v1.0   27/04/2018 Changed version number to match MATLAB Central
% v0.2.3 16/03/2018 Changed SDF_CHANNEL_HDR.overloaded template entry to
% deal with error in an overloaded file. The value '3' (out of 0-1 range)
% was presen for a file which exhibited 'overloaded' on the instrument
% screen.
% v0.2.2 06/04/2017 Changed the untested decimated time data field name to decimatedRealValue and
%   decimated ImaginaryValue, to reduce confusion associated with using minimum and maximum Real and 
%   Imaginary values (see B-37 of standard)
% v0.2.1 04/04/2017 Reformatted under Matlab, and recast everything to
%   double to avoid values being quantised.
% v0.2 - 03/04/2017 Corrected major Matlab/Octave syntax issues:
%            - Use of '_' at start of variable names replaced with 'zz_'
%            - 'data._SDF_ADDTL' renamed 'data.SDF_ADDTL'
%            - Replaced '!' with '~' for the 'not' statements
%            - Replaced double quotes (") with single quotes (') as appropriate
% v0.1 - 31/03/2017 Initial alpha release

%   This file is part of SDF Importer.
% 
% Copyright 2016-2018 Commonwealth of Australia, represented by the Department of Defence.
%
% Redistribution and use in source and binary forms, with or without modification, are 
% permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this list of 
%    conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice, this list 
%    of conditions and the following disclaimer in the documentation and/or other materials 
%    provided with the distribution.
%
% 3. Neither the name of the copyright holder nor the names of its contributors may be used 
%    to endorse or promote products derived from this software without specific prior 
%    written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY 
% AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% Initialise Variables
data.ext.traceList = [{},{}];
data.SDF_ADDTL.DATA_HDR_sizes = 0; % Array of recordSize entries from the DATA_HDR header
data.SDF_ADDTL.VECTOR_HDR_sizes = 0; % Array of recordSize entries from the VECTOR_HDR header
data.SDF_ADDTL.CHANNEL_HDR_sizes = 0; % Array of recordSize entries from the CHANNEL_HDR header
data.SDF_ADDTL.UNIQUE_sizes = 0; % Array of recordSize entries from the UNIQUE header
data.SDF_ADDTL.SCAN_STRUCT_sizes = 0; % Array of recordSize entries from the SCAN_STRUCT header
data.SDF_ADDTL.SCAN_BIG_sizes = 0; % Array of recordSize entries from the SCAN_BIG header
data.SDF_ADDTL.SCANS_VAR_sizes = 0; % Array of recordSize entries from the SANS_VAR header
data.SDF_ADDTL.COMMENT_sizes = 0;  % Array of recordSize entries from the COMMENT header
data.SDF_ADDTL.filename=filename; % Store filename for reference in error messages and warnings

[fid,message] = fopen(filename,'r'); % Open file in read only mode to avoid corruption
if fid == -1
    disp([ mfilename ': ERROR opening file ' filename ])
    error(message)
end
file_raw = fread(fid); % Read file marker to make sure it is SDF file
fclose(fid); % Close file

%% Confirm the file contains an SDF header by looking for 'B\0' as first two bytes of the file.
if strcmp(char(file_raw(1:2)),[ 'B'; char(0) ])
    
    
    % Import SDF_FILE_HDR
    data.SDF_FILE_HDR = SDF_Interpret_Template('SDF_FILE_HDR', file_raw, 3);
    
    % Import SDF_MEAS_HDR
    data.SDF_MEAS_HDR = SDF_Interpret_Template('SDF_MEAS_HDR', file_raw,data.SDF_FILE_HDR.recordSize+3);
    
    % Estimate the version of the SDF file
    % Based on the size of MEAS_HDR, we deduce the SDF Format version
    switch data.SDF_MEAS_HDR.recordSize
        case 102
            data.SDF_ADDTL.SDF_version = 1.0;
        case 140
            data.SDF_ADDTL.SDF_version = 2.0;
        case 156
            data.SDF_ADDTL.SDF_version = 3.0;
        otherwise % Catch-all case incase the standard is modified
            error('Unrecognised SDF version: data.SDF_MEAS_HDR.recordSize = %d',...
                data.SDF_MEAS_HDR.recordSize);
    end
    
    % Import SDF_DATA_HDR
    if (data.SDF_FILE_HDR.num_of_DATA_HDR_record>0)
        for n = 1:data.SDF_FILE_HDR.num_of_DATA_HDR_record % Iterated over all the DATA_HDR records
            % Import the nth 'DATA' record
            data.SDF_DATA_HDR(n) = ...
                SDF_Interpret_Template('SDF_DATA_HDR', ...
                file_raw,data.SDF_FILE_HDR.offset_of_DATA_HDR_record + 1 + ...
                sum(data.SDF_ADDTL.DATA_HDR_sizes));
            % Update vector of record sizes
            data.SDF_ADDTL.DATA_HDR_sizes(n) = data.SDF_DATA_HDR(n).recordSize;
        end
    end
    
    % Import SDF_VECTOR_HDR
    if (data.SDF_FILE_HDR.num_of_VECTOR_record>0)
        for n = 1:data.SDF_FILE_HDR.num_of_VECTOR_record
            % Import the nth 'VECTOR'record
            data.SDF_VECTOR_HDR(n) = ...
                SDF_Interpret_Template('SDF_VECTOR_HDR', ...
                file_raw,data.SDF_FILE_HDR.offset_of_VECTOR_record + 1 + ...
                sum(data.SDF_ADDTL.VECTOR_HDR_sizes));
            % Update vector of record sizes
            data.SDF_ADDTL.VECTOR_HDR_sizes(n) = data.SDF_VECTOR_HDR(n).recordSize;
        end
    end
    
    % Import SDF_CHANNEL_HDR
    if (data.SDF_FILE_HDR.num_of_CHANNEL_record>0)
        for n = 1:data.SDF_FILE_HDR.num_of_CHANNEL_record
            % Import the nth 'CHANNEL'record
            data.SDF_CHANNEL_HDR(n) = ...
                SDF_Interpret_Template('SDF_CHANNEL_HDR', ...
                file_raw,data.SDF_FILE_HDR.offset_of_CHANNEL_record + 1 + ...
                sum(data.SDF_ADDTL.CHANNEL_HDR_sizes));
            % Update vector of record sizes
            data.SDF_ADDTL.CHANNEL_HDR_sizes(n) = data.SDF_CHANNEL_HDR(n).recordSize;
        end
    end
    
    % Import UNIQUE records
    if (data.SDF_FILE_HDR.num_of_UNIQUE_record > 0)
        for n = 1:data.SDF_FILE_HDR.num_of_UNIQUE_record
            % Import the nth 'UNIQUE'record
            data.SDF_UNIQUE(n) = ...
                SDF_Interpret_Template('SDF_UNIQUE', ...
                file_raw,data.SDF_FILE_HDR.offset_of_UNIQUE_record + 1 + ...
                sum(data.SDF_ADDTL.UNIQUE_sizes));
            % Update vector of record sizes
            data.SDF_ADDTL.UNIQUE_sizes(n) = data.SDF_UNIQUE(n).recordSize;
        end
    end
    
    % Import SCAN_STRUCT_record
    if (data.SDF_FILE_HDR.num_of_SCAN_STRUCT_record>0)
        for n = 1:data.SDF_FILE_HDR.num_of_SCAN_STRUCT_record
            % Import the nth 'SCAN_STRUCT'record
            data.SDF_SCAN_STRUCT(n) = ...
                SDF_Interpret_Template('SDF_SCAN_STRUCT', ...
                file_raw,data.SDF_FILE_HDR.offset_of_SCAN_STRUCT_record + 1 + ...
                sum(data.SDF_ADDTL.SCAN_STRUCT_sizes));
            % Update vector of record sizes
            data.SDF_ADDTL.SCAN_STRUCT_sizes(n) = data.SDF_SCAN_STRUCT(n).recordSize;
            % Add in interpretation of scan struct raw_data for HP89410 and HP9440 (Applies to all?)
            y = reshape(data.SDF_SCAN_STRUCT.zz_raw_data,data.SDF_SCAN_STRUCT.num_of_scan,...
                length(data.SDF_SCAN_STRUCT.zz_raw_data)/data.SDF_SCAN_STRUCT.num_of_scan);
            eval(['data.SDF_SCAN_STRUCT.scanVar(n,:) = SDF_' data.SDF_SCAN_STRUCT.zz_scanVar_type ...
                '_vector(y);']);
            
        end
    end
    
    % Check the file version support SCAN_BIG, SCAN_VAR, and COMMENTS ( >= 3.0)
    if(data.SDF_ADDTL.SDF_version  >= 3.0)
        % Import all SCAN_BIG
        % Check there are SCANG_BIG records
        if (data.SDF_FILE_HDR.num_of_SCAN_BIG_RECORD>0)
            warning([data.SDF_ADDTL.filename ': SCAN_BIG_RECORD, entering untested section of code.']);
            % Iterate over the SCAN_BIG records
            for n = 1:data.SDF_FILE_HDR.num_of_SCAN_BIG_RECORD
                % Import the nth 'SCAN_BIG'record
                data.SDF_SCAN_BIG(n) = ...
                    SDF_Interpret_Template('SDF_SCAN_BIG', ...
                    file_raw,data.SDF_FILE_HDR.offset_of_SCAN_BIG_record + 1 + ...
                    sum(data.SDF_ADDTL.SCAN_BIG_sizes));
                % Update vector of record sizes
                data.SDF_ADDTL.SCAN_BIG_sizes(n) = data.SDF_SCAN_BIG(n).recordSize;
            end
            
            
            % We have not implemented the SDF_SCANS_VAR importation.
            % We make the assumption that the SCANS_VAR is located just after the SCAN_BIG
            % section, so we will define an additonal variable:
            % SDF_ADDTL.offset_of_SCAN_VAR_record
            warning([data.SDF_ADDTL.filename ': SDF_SCANS_VAR importation has not been implemented.']);
            warning('data.SDF_ADDTL.offset_of_SCANS_VAR_record has been created');
            % Initialise to zero in case of empty set.
            data.SDF_ADDTL.SCANS_VAR_sizes = 0;
            % Calculate offset based on previous dataset
            data.SDF_ADDTL.offset_of_SCANS_VAR_record = ...
                data.SDF_FILE_HDR.offset_of_SCAN_BIG_record + ...
                sum(data.SDF_ADDTL.SCAN_BIG_sizes);
            
            if (data.SDF_FILE_HDR.num_of_COMMENT_record>0)
                warning([data.SDF_ADDTL.filename ': SDF_COMMENT_HDR importation has not been implemented.']);
                warning('data.SDF_ADDTL.offset_of_COMMEENT_record has been created');
                %% We have not implemented the SDF_COMMENT_HDR importation.
                % Initialise to zero in case of empty set.
                data.SDF_ADDTL.COMMENT_sizes = 0;
                % Calculate offset based on previous dataset
                data.SDF_ADDTL.offset_of_COMMENT_record = ...
                    data.SDF_FILE_HDR.offset_of_SCANS_VAR_record + ...
                    sum(data.SDF_ADDTL.SCANS_VAR_sizes);
            end %num_of_COMMENT_record
        end%num_of_SCAN_BIG_RECORD
    end
    
    % Import all YDATA
    data = SDF_YDATA_Process(data,file_raw);
    
    % Import all XDATA if necessary
    data = SDF_XDATA_Process(data,file_raw);
else %if strcmp(char(file_raw(1:2)),[ 'B'; char(0) ])
    error(['The input file: ' filename 'is not an SDF file']);
end  %if strcmp(char(file_raw(1:2)),[ 'B'; char(0) ]), else
end % function [ data file_raw] = SDF_import( filename )

function [retval] = SDF_Interpret_Template(template, file_raw, start_address)
%SDF_INTERPRET_TEMPLATE - Iterates over the elements in the array of structures
% referred to by the string in template, reading input byte array 'file_raw',
% from index 'start_address', returning the new structure in 'retval'
%
%  This code is an implementation of the SDF standard as defined in:
%  'User's Guide - Standard Data Format Utilities, Version B.02.01'
%  by Agilent Technologies
%  Manufacturing Part Number: 5963-1715
%  Printed in USA
%  December 1994
%  © Copyright 1989, 1991-94, 2005 Agilent Technologies, Inc.
%
% Syntax: [retval] = SDF_Interpret_Template(template, file_raw, start_address)
%
% Inputs:
%   template - String of the desired header, options are:
%                 - SDF_FILE_HDR:    Provides an index to the file.
%                 - SDF_MEAS_HDR:    Contains settings of measurement parameters.
%                 - SDF_DATA_HDR:    Tells you how reconstruct one block of measurement results
%                                    (x- and y-axis values for every point of every trace).
%                 - SDF_VECTOR_HDR:  Tells you which channel (or pair of channels) provided data
%                                    for a single trace.
%                 - SDF_CHANNEL_HDR: Contains channel-specific information for one channel
%                                    used in the measurement.
%                 - SDF_SCAN_STRUCT: Tells you how vectors are organized in the Y-axis Data
%                                    record when the measurement includes multiple scans of data.
%                 - SDF_SCANS_BIG:   Extended scan header which tells how vectors are organized
%                                    in the Y-axis data record when the measurement may include
%                                    more than 32767 scans of data.
%                 - SDF_SCANS_VAR:   Contains a number that identifies every scan
%                                    (scan time, RPM, or scan number).
%                 - SDF_COMMENT_HDR: Contains text information associated with the file.
%                 - SDF_XDATA_HDR:   Contains the x-axis data needed to reconstruct any trace.
%                 - SDF_YDATA_HDR:   Contains the y-axis data needed to reconstruct any trace.
%                 - SDF_UNIQUE:      Makes the SDF file format flexible. The eight common record
%                                    types define parameters that are common to many instruments
%                                    and systems. However, a particular instrument or system may
%                                    need to save and recall the states of additional parameters.
%                                    These states can reside in Unique records.
%                 - SDF_UNIT:        Contains eng. units & scaling information for the traces.
%                 - SDF_WINDOW:      Contains windowing information for frequency domain traces.
%   file_raw - Raw byte array of the SDF file which is being intepreted.
%   start_address - array index of the header's first element.
%
% Outputs:
%   retval - Header structure with data populated as specified in the template.
%
% Other m-files required: SDF_template.m
% Subfunctions: SDF_DataType, SDF_uintX, SDF_int8, SDF_char, SDF_short, SDF_long, SDF_float,
%   SDF_double, SDF_struct,
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

% Get the retval Template
tmplt = SDF_template(template);

% Use start address as a flag to identify if file_raw is all the imported file,
% or just the section of data relevant to SDF_UNIT or SDF_WINDOW templates
if (start_address >0) % The case for HDR template
    % Interpret the header length bytes.
    recordSize = SDF_long(file_raw(start_address+2:start_address+2+3));
    % Get only the data relevant to the HDR being imported
    data = file_raw(start_address:recordSize+start_address-1);
else % The case for  SDF_UNIT or SDF_WINDOW templates
    % The relevant data will have already been selected before being passed.
    data = file_raw;
    recordSize = length(file_raw); % Allows all data to be processed.
end

% idx is used to track the largest 'Binary Index' referred to by the header.
% It is necessary because to maintain the same 'Field Index' as the SDF doc.
% Additional 'Field Index' values were created to account for struct arrays.
idx = 0;
retval = []; % Initiate retval (return value) variable
% Iterate over the template entries to extract data and place in retval
for cnt = tmplt
    if cnt.BinaryIndex(end)  <= recordSize % Confirm record has not ended
        retval = SDF_DataType(cnt,data,retval); % Convert single element of data (not vector)
        idx = max([idx cnt.BinaryIndex(end)+1]); % Track the byte after the current element
    end
end

% Iterate over entries reproduce 'Old' variables as 'new' if they don't exist
for cnt = tmplt
    if strcmp(cnt.FieldName(end-2:end),'Old') % Identify 'Old' variables from the template
        if (~isfield(retval, cnt.FieldName(1:end-3))) % Identify non-existant variables
            % Copy the old address to the new one variable location
            eval(['retval.' cnt.FieldName(1:end-3) ' = retval.' cnt.FieldName ';']);
        end
    end
end

for cnt = tmplt % Iterate over entries to populate any entries with a table associated.
    if cnt.BinaryIndex(end)  <= recordSize
        if (size(cnt.Table,1) ~= 0) % Populate Fields with Table associated
            eval(['retval.zz_' cnt.FieldName ' = cnt.Table{' ...
                'find(strcmp(cnt.Table,num2str(retval.' cnt.FieldName '))),2};']);
        end
    end
end

% If there is data remaining in header, dump into _raw_data
if (start_address  ~= 0)
    tmp = data(idx:retval.recordSize);
    if (numel(tmp) ~= 0) % Only create variable if there is data
        retval.zz_raw_data = tmp;
    end
end
% End of the function
end

function [retval] = SDF_YDATA_Process(data,file_raw)
%SDF_YDATA_PROCESS - Identifies how y-data is stored in file_raw, and then extracts it to a
%  standardised structure.
%
%  This code is an implementation of the SDF standard as defined in:
%  'User's Guide - Standard Data Format Utilities, Version B.02.01'
%  by Agilent Technologies
%  Manufacturing Part Number: 5963-1715
%  Printed in USA
%  December 1994
%  © Copyright 1989, 1991-94, 2005 Agilent Technologies, Inc.
%
% Syntax: [retval] = SDF_YDATA_Process(data,file_raw)
%
% Inputs:
%   data - Structure containing the populated SDF header information.
%   file_raw - Raw byte array of the SDF file which is being intepreted.
%
% Outputs:
%   retval - Input variable 'data' with the y-data added
%
% Other m-files required: none
% Subfunctions: capitaliseFirstLetters, dataTypeByteSize, SDF_DataType,
%   SDF_uintX, SDF_int8, SDF_char, SDF_short, SDF_long, SDF_float, SDF_double, SDF_struct,
%   SDF_uintX_vector, SDF_short_vector, SDF_long_vector, SDF_float_vector, SDF_double_vector,
%   SDF_DATA_SizesY, SDF_extract_labels, SDF_Interpret_Template, SDF_Multi_ChannelY,
%   SDF_Multi_TraceX, SDF_Multi_TraceY, SDF_Single_TraceY, SDF_YDATA_Process, SDF_Y_scale.
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

% Calculate the offsets for each new set of Y-data
data_offsets = [0 data.SDF_ADDTL.DATA_HDR_sizes(1:end-1)];

% Import the YDATA header and its raw_data based on the DATA_HDR offsets.
data.SDF_YDATA_HDR = ...
    SDF_Interpret_Template ('SDF_YDATA_HDR',file_raw,data.SDF_FILE_HDR.offset_of_YDATA_record+1 );

% Find the size of each vector/data block so we can work at extracting the correct vectors.
data = SDF_DATA_SizesY(data);

% Based on the scan structure, interpret data differently.
if data.SDF_FILE_HDR.num_of_SCAN_STRUCT_record % If there is a scan structure check what type
    switch data.SDF_SCAN_STRUCT.zz_scan_type
        case 'scan'
            y_idx = 0; % Initialise the index counter for the block of data
            % Scan #, Data #, Vector #,
            for n_s = 1 : data.SDF_SCAN_STRUCT.num_of_scan    %Iterate over all the scans
                for n_d = 1 : sum(data.SDF_ADDTL.is_scan_data) %Assume sequential, iterate datasets
                    [data,y_idx] = SDF_Multi_TraceY(data,n_s,n_d,y_idx);
                end %n_
            end %n_s
            % Assume the non-scan data is after the scan structured data, import
            % the non-scanned data starting with the first non-scan datatype
            % Note we assume the scans are all before the non-scan.
            for n_ns = sum(data.SDF_ADDTL.is_scan_data)+1:length(data.SDF_ADDTL.is_scan_data)
                [data,y_idx] = SDF_Multi_ChannelY(data,n_ns,y_idx); % For some reason 'data' was 'tmp'
            end %n_ns
            
        case 'depth'
            warning([data.SDF_ADDTL.filename ': An untested section of code was executed in '...
                'SDF_YDATA_Process. data.SDF_SCAN_STRUCT.zz_scan_type = ''depth''']);
            % Data #, Scan #, Vector #
            y_idx = 0; % Initialise the index counter for the block of data
            % Scan #, Data #, Vector #,
            for n_d = 1 : sum(data.SDF_ADDTL.is_scan_data)  %Iterate over all the scans
                for n_s = 1 : data.SDF_SCAN_STRUCT.num_of_scan %Assume sequential, iterate datasets
                    [data,y_idx] = SDF_Multi_TraceY(data,n_s,n_d,y_idx);
                end %n_
            end %n_s
            % Assume the non-scan data is after the scan structured data, import
            % the non-scanned data starting with the first non-scan datatype
            % Note we assume the scans are all before the non-scan.
            for n_ns = sum(data.SDF_ADDTL.is_scan_data)+1:length(data.SDF_ADDTL.is_scan_data)
                [data,y_idx] = SDF_Multi_ChannelY(data,n_ns,y_idx); % For some reason 'data' was 'tmp'
            end %n_ns
        otherwise
            % Do nothing
            warning([data.SDF_ADDTL.filename ': An unreachable section of code was executed in '...
                'SDF_YDATA_Process. data.SDF_SCAN_STRUCT.zz_scan_type  ~= ''depth'' | ''scan''']);
    end % switch data.SDF_SCAN_STRUCT.zz_scan_type
else % There is no scan struct, read y trace
    data = SDF_Single_TraceY(data); % Old way for reading a single trace
end % if data.SDF_FILE_HDR.num_of_SCAN_STRUCT_record
% Return data
retval = data;
% End of function
end




function [retval] = SDF_DataType(cnt,raw_data,data)
%SDF_DATATYPE - Interprets a sequence of bytes within data, return as datatype
%  cnt.DataType. retval is updated with the interpreted value.
%
%  This code is an implementation of the SDF standard as defined in:
%  'User's Guide - Standard Data Format Utilities, Version B.02.01'
%  by Agilent Technologies
%  Manufacturing Part Number: 5963-1715
%  Printed in USA
%  December 1994
%  © Copyright 1989, 1991-94, 2005 Agilent Technologies, Inc.
%
% Syntax: [retval] = SDF_DataType(cnt,raw_data,data)
%
% Inputs:
%   cnt - Nth sub-structure from SDF template structure.
%   raw_data - Raw byte array of the SDF file sub-section which is being intepreted.
%   data - Structure containing the populated SDF header information.
%
% Outputs:
%   retval - Input variable 'data' with the interpreted data added
%
% Other m-files required: none
% Subfunctions: SDF_uintX, SDF_int8, SDF_char, SDF_short, SDF_long, SDF_float, SDF_double,
%               SDF_struct.
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

% Preload the existing data into the new variable to be returned.
retval = data;
% Catch the char and struct fields by looking at their first four characters
switch cnt.DataType(1:4)
    % Process all strings. Single strings are 'int8', and dealt with by the
    % SDF_char function. If the standard is updated to include single chars
    % which are displayed as chars, a different implementation is needed.
    case {'char'}
        eval(['retval.' cnt.FieldName ' = SDF_char' ...
            '(raw_data(cnt.BinaryIndex(1):cnt.BinaryIndex(end)));']);
        % Process the structures in their own way
    case {'stru'} %struct
        eval(['retval.' cnt.FieldName ' = SDF_struct(''' cnt.RangeUnits(5:end)...
            ''',raw_data(cnt.BinaryIndex(1):cnt.BinaryIndex(end)));']);
        % Process all standard datatypes as per their name
    otherwise
        eval(['retval.' cnt.FieldName ' = SDF_' cnt.DataType ...
            '(raw_data(cnt.BinaryIndex(1):cnt.BinaryIndex(end)));']);
end
% End of function
end


function [retval] = SDF_uintX (arg)
%SDF_UINTX - converts a 1D array of N = 2:2:8 bytes to an unsigned integer

%  concatenanting the bits in reverse byte order. (To account for SDF byte order)
%
% Syntax: [retval] = SDF_uintX( arg )
%
% Inputs:
%   arg - 1D array of bytes of length 2 to 8.
%
% Outputs:
%   retval - Unsigned integer of N bits, where N = 8*length(arg);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

% Reverse bytes
arg = arg(end:-1:1)';
% Get number of bytes (2:2:8)
bytes = length(arg);
eval(['retval = uint' num2str(bytes * 8) '(0);']);
% Create index of bytes
n = (1:bytes)';
% Use vector operations to scale the bytes appropriately
try
  retval = retval+arg * (2.^((n-1) * 8));
catch
  retval = uint32(0);
  retval = retval+arg * (2.^((n-1) * 8));
end
% End function
end

function [retval] = SDF_uintX_vector (arg)
%SDF_UINTX_VECTOR - converts a 2D array of N x M where N = 2:2:8 bytes to an unsigned integer
%  concatenanting the bits in reverse byte order. (To account for SDF byte order)
%
% Syntax: [retval] = SDF_uintX_vector( arg )
%
% Inputs:
%   arg - 2D array of bytes of 2 to 8 rows, by M columns.
%
% Outputs:
%   retval - Array of length M of unsigned integer of N bits, where N = 8*size(arg,1);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

% Reverse bytes
arg = arg(end:-1:1,:);
% Get number of bytes (2:2:8)
bytes = size(arg,1);
col = size(arg,2);
eval(['retval = uint' num2str(bytes * 8) '(0);']);
retval = repmat(retval,1,col);
eval(['arg = uint' num2str(bytes * 8) '(arg);']);
% Create index of bytes

eval(['n = uint' num2str(bytes * 8) '(((1:bytes)-1) * 8);']);
n = repmat(transpose(n),1,col);
% Use vector operations to scale the bytes appropriately
retval = (2.^n).* arg; % Scale the bytes by powers of two
retval = eval(['uint' num2str(bytes * 8) '(sum(retval,1))']); % Sum columns to get the vector.
% End function
end

function [output] = SDF_char(arg)
%SDF_CHAR - converts an array of N bytes to a string in
%  reverse byte order. (To account for SDF byte order)
%
% Syntax: [retval] = SDF_char( arg )
%
% Inputs:
%   arg - 2D array of bytes of length N
%
% Outputs:
%   retval - Array of char of length N;
%
% Other m-files required: none
% Subfunctions: SDF_int8
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

% Reverse bytes
arg = arg';
% Check if sting or single char (int8?)
if length(arg)>1
    % Convert to a text string
    output = char(arg);
    % Remove unprintable chars
    output(~ismember(double(output),32:255)) = '';
else
    % Convert the char to an int8
    output = SDF_int8(arg);
end
% End function
end

function [retval] = SDF_int8(arg)
%SDF_INT8 - converts a 1D array of N = 1 bytes to a signed 8-bit integer.
%
% Syntax: [retval] = SDF_int8( arg )
%
% Inputs:
%   arg - 1D array of bytes of length 1.
%
% Outputs:
%   retval - arg interpreted as an int8;
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

if  bitget(arg,8)%Convert if 8th bit is 1 (negative number)
    arg = arg - 254;
end
retval = arg;
end

function [retval] = SDF_short(arg)
%SDF_SHORT - converts a 1D array of N = 2 bytes to a short (signed 16 bit integer)
%  concatenanting the bits in reverse byte order. (To account for SDF byte order)
%
% Syntax: [retval] = SDF_short( arg )
%
% Inputs:
%   arg - 1D array of bytes of length 2.
%
% Outputs:
%   retval - arg interpreted as an int16;
%
% Other m-files required: none
% Subfunctions: SDF_uintX
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

% Merge send the bytes to be merged in reverse order and typecast them.
retval = double(int32(typecast(SDF_uintX(arg),'int16')));
end

function [output] = SDF_long(arg)
%SDF_LONG - converts a 1D array of N = 4 bytes to a short (signed 32 bit integer)
%  concatenanting the bits in reverse byte order. (To account for SDF byte order)
%
% Syntax: [retval] = SDF_long( arg )
%
% Inputs:
%   arg - 1D array of bytes of length 4.
%
% Outputs:
%   retval - arg interpreted as an int32;
%
% Other m-files required: none
% Subfunctions: SDF_uintX
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

% Merge send the bytes to be merged in reverse order and typecast them.
output = double(typecast(SDF_uintX(arg),'int32'));
end

function [retval] = SDF_float(arg)
%SDF_FLOAT - converts a 1D array of N = 4 bytes to a float (signed 32 bit floating point precision,
%  also known as a single) concatenanting the bits in reverse byte order.
%  (To account for SDF byte order)
%
% Syntax: [retval] = SDF_float( arg )
%
% Inputs:
%   arg - 1D array of bytes of length 4.
%
% Outputs:
%   retval - arg interpreted as an float (single);
%
% Other m-files required: none
% Subfunctions: SDF_uintX
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

% Merge send the bytes to be merged in reverse order and typecast them.
retval = double(typecast(SDF_uintX(arg),'single'));
end

function [retval] = SDF_double(arg)
%SDF_DOUBLE - converts a 1D array of N = 8 bytes to a double (signed 64 bit floating point
%  precision) concatenanting the bits in reverse byte order. (To account for SDF byte order)
%
% Syntax: [retval] = SDF_double( arg )
%
% Inputs:
%   arg - 1D array of bytes of length 8.
%
% Outputs:
%   retval - arg interpreted as a double;
%
% Other m-files required: none
% Subfunctions: SDF_uintX
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

% Merge send the bytes to be merged in reverse order and typecast them.
retval = cast(SDF_uintX(arg),'double');
end


function [retval] = SDF_short_vector(arg)
%SDF_SHORT_VECTOR - converts a 2D array of (N = 2) x M bytes to an array of short
%  (signed 16 bit integers) concatenanting the bits in reverse byte order.
%  (To account for SDF byte order)
%
% Syntax: [retval] = SDF_short_vector( arg )
%
% Inputs:
%   arg - 2D array of bytes of length 2xM.
%
% Outputs:
%   retval - Array of length M of shorts (integer of 16 bits);
%
% Other m-files required: none
% Subfunctions: SDF_uintX_vector
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

% Merge send the bytes to be merged in reverse order and typecast them.
retval = double(int32(typecast(SDF_uintX_vector(arg),'int16')));
end

function [output] = SDF_long_vector(arg)
%SDF_LONG_VECTOR - converts a 2D array of (N = 4) x M bytes to an array of short
%  (signed 16 bit integers) concatenanting the bits in reverse byte order.
%  (To account for SDF byte order)
%
% Syntax: [retval] = SDF_long_vector( arg )
%
% Inputs:
%   arg - 2D array of bytes of length 4xM.
%
% Outputs:
%   retval - Array of length M of shorts (integer of 32 bits);
%
% Other m-files required: none
% Subfunctions: SDF_uintX_vector
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

% Merge send the bytes to be merged in reverse order and typecast them.
output = double(typecast(SDF_uintX_vector(arg),'int32'));
end

function [retval] = SDF_float_vector(arg)
%SDF_FLOAT_VECTOR - converts a 2D array of (N = 4) x M bytes to an array of short
%  (signed 16 bit integers) concatenanting the bits in reverse byte order.
%  (To account for SDF byte order)
%
% Syntax: [retval] = SDF_float_vector( arg )
%
% Inputs:
%   arg - 2D array of bytes of length 4xM.
%
% Outputs:
%   retval - Array of length M of float (single - 32 bits of floating pint precision);
%
% Other m-files required: none
% Subfunctions: SDF_uintX_vector
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

% Merge send the bytes to be merged in reverse order and typecast them.
retval = double(typecast(SDF_uintX_vector(arg),'single'));
end

function [retval] = SDF_double_vector(arg)
%SDF_DOUBLE_VECTOR - converts a 2D array of (N = 8) x M bytes to an array of short
%  (signed 16 bit integers) concatenanting the bits in reverse byte order.
%  (To account for SDF byte order)
%
% Syntax: [retval] = SDF_DOUBLE_vector( arg )
%
% Inputs:
%   arg - 2D array of bytes of length 8xM.
%
% Outputs:
%   retval - Array of length M of double (double - 64 bits of floating pint precision);
%
% Other m-files required: none
% Subfunctions: SDF_uintX_vector
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

% Merge send the bytes to be merged in reverse order and typecast them.
retval = typecast(SDF_uintX_vector(arg),'double');
end

function [retval] = SDF_struct(struct_template,arg)
%SDF_STRUCT - converts a 1D array of bytes to a SDF_UNIT or SDF_WINDOW structure.
%
% Syntax: [retval] = SDF_struct( arg )
%
% Inputs:
%   struct_template - A text string identifying the relectan structure to apply, options:
%                       'SDF_UNIT'
%                       'SDF_WINDOW'
%   arg - 1D array of bytes.
%
% Outputs:
%   retval - arg interpreted SDF_UNIT of SDF_WINDOW(single);
%
% Other m-files required: none
% Subfunctions: SDF_Interpret_Template
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

% Check which struct needs to be processed, return error if not implemented
switch struct_template
    case {'SDF_UNIT'}
        retval = SDF_Interpret_Template('SDF_UNIT', arg, 0);
    case {'SDF_WINDOW'}
        retval = SDF_Interpret_Template('SDF_WINDOW', arg, 0);
    otherwise
        error('This is a struct which is not implemented: %s',struct_template)
end
end


function [retval] = SDF_Single_TraceY(data)
%SDF_SINGLE_TRACEY - Reconstructs the Y-data for the simplest type of SDF file, which contains a
%  single trace.
%  This code is an implementation of the SDF standard as defined in:
%  'User's Guide - Standard Data Format Utilities, Version B.02.01'
%  by Agilent Technologies
%  Manufacturing Part Number: 5963-1715
%  Printed in USA
%  December 1994
%  © Copyright 1989, 1991-94, 2005 Agilent Technologies, Inc.
%
% Syntax: [retval] = SDF_Single_TraceY(data)
%
% Inputs:
%   data - Extracted SDF headers, containing the raw y-data.
%
% Outputs:
%   retval - Extracted SDF headers, and processed y-axis traces.
%
% Other m-files required: none
% Subfunctions: capitaliseFirstLetters, dataTypeByteSize, SDF_DataType,
%   SDF_uintX_vector, SDF_short_vector, SDF_long_vector, SDF_float_vector, SDF_double_vector,
%   SDF_DATA_SizesY, SDF_extract_labels, SDF_Y_scale.
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

% Perform the import of data for a single trace.
var_name = 'data.ext.s0d1v1.y'; % Define variable name
% var_size is the size (in bytes) of the Y-axis Data record’s data block
for n = 1:length(data.SDF_DATA_HDR)
    % If there is a single trace immediately sort data chunks are individual datapoint values.
    data.SDF_YDATA_HDR.zz_raw_data = ...
        reshape(data.SDF_YDATA_HDR.zz_raw_data, ...
        data.SDF_ADDTL.YDATA_var_size,...
        length(data.SDF_YDATA_HDR.zz_raw_data)/...
        data.SDF_ADDTL.YDATA_var_size);
    % Merge the bytes within the rows to allow byte conversion
    eval(['data.SDF_YDATA_HDR.zz_data = '...
        'SDF_' data.SDF_DATA_HDR(n).zz_ydata_type '_vector(data.SDF_YDATA_HDR.zz_raw_data);']);
    % Abbreviate some variable names for neatness
    dataHdr = data.SDF_DATA_HDR;
    vctHdr = data.SDF_VECTOR_HDR(dataHdr.first_VECTOR_recordNum+1);
    chHdr = data.SDF_CHANNEL_HDR;
    tmp = data.SDF_YDATA_HDR.zz_data;
    % Scale the converted data chuncks from volts to engineering units, and spit by yPerPoint and
    % by isComplex
    eval([var_name ' = SDF_Y_scale(chHdr,dataHdr,vctHdr,tmp);']); % scale the results
    % Extract user-readable data assiciated with the processed channels/axes.
    data = SDF_extract_labels(data,var_name);
    retval = data;
end
end

function [retval,y_idx] = SDF_Multi_TraceY(data,n_s,n_d,y_idx)
%SDF_MULTI_TRACEY - Reconstructs the scan Y-data for the datasets which are part of a multiscan file.
%  This code is an implementation of the SDF standard as defined in:
%  'User's Guide - Standard Data Format Utilities, Version B.02.01'
%  by Agilent Technologies
%  Manufacturing Part Number: 5963-1715
%  Printed in USA
%  December 1994
%  © Copyright 1989, 1991-94, 2005 Agilent Technologies, Inc.
%
% Syntax: [retval] = SDF_Multi_TraceY(data)
%
% Inputs:
%   data - Extracted SDF headers, containing the raw y-data.
%   n_s - Scan number (must be greater than zero)
%   n_d - Data header number (must be greater than zero)
%   y_idx - index of the first byte relevant to the Multi Trace input data.SDF_YDATA_HDR.zz_raw_data
%
% Outputs:
%   retval - Extracted SDF headers, and processed y-axis traces.
%
% Other m-files required: none
% Subfunctions: capitaliseFirstLetters, dataTypeByteSize, SDF_DataType,
%   SDF_uintX_vector, SDF_short_vector, SDF_long_vector, SDF_float_vector, SDF_double_vector,
%   SDF_DATA_SizesY, SDF_extract_labels, SDF_Y_scale.
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

%Calculate the number of vectors for this datatype
num_v = data.SDF_DATA_HDR(n_d).total_rows * ...
    data.SDF_DATA_HDR(n_d).total_cols;
% Read the block size we have calculated
blk_size = data.SDF_ADDTL.YDATA_block_sizes(n_d); %Includes yPerPoint and isComplex
% Import all the vectors reshaped to be in rows
for n_v = 1:num_v
    strt = y_idx+(n_v-1) * blk_size+1; % Start index of data of interest
    stp = y_idx+(n_v) * blk_size; % Stop index of data of interest
    
    
    % Define the variable name
    var_name = ['data.ext.s' num2str(n_s) ...
        'd' num2str(n_d) 'v' num2str(n_v) '.y'];
    % Abbreviate some variable names for neatness
    dataHdr = data.SDF_DATA_HDR(n_d);
    vctHdr = data.SDF_VECTOR_HDR(dataHdr.first_VECTOR_recordNum+n_v);
    chHdr = data.SDF_CHANNEL_HDR;
    % Place the relevant data into reshaped tmp array vector
    tmp = reshape(data.SDF_YDATA_HDR.zz_raw_data(strt:stp),...
        data.SDF_ADDTL.YDATA_var_size(n_d), ...
        blk_size/int32(data.SDF_ADDTL.YDATA_var_size(n_d)));
    % Convert the vector to datatype and scale in single command.
    eval([var_name ' = SDF_Y_scale(chHdr,dataHdr,vctHdr,' ...
        'SDF_' data.SDF_DATA_HDR(n_d).zz_ydata_type '_vector(tmp));']); % scale the results
    % Extract user-readable data assiciated with the processed channels/axes.
    data = SDF_extract_labels(data,var_name);
end
% Increment the index counter
y_idx = stp;
% Return data
retval = data;
end




function [retval,y_idx] = SDF_Multi_ChannelY(data,n_ns,y_idx)
%SDF_MULTI_CHANNELY - Reconstructs non-scan Y-data for the datasets which are part of a multiscan file.
%  This code is an implementation of the SDF standard as defined in:
%  'User's Guide - Standard Data Format Utilities, Version B.02.01'
%  by Agilent Technologies
%  Manufacturing Part Number: 5963-1715
%  Printed in USA
%  December 1994
%  © Copyright 1989, 1991-94, 2005 Agilent Technologies, Inc.
%
% Syntax: [retval] = SDF_Multi_ChannelY(data)
%
% Inputs:
%   data - Extracted SDF headers, containing the raw y-data.
%   n_ns - Data header number (must be greater than zero)
%   y_idx - index of the first byte relevant to the Multi Channel input data.SDF_YDATA_HDR.zz_raw_data
%
% Outputs:
%   retval - Extracted SDF headers, and processed y-axis traces.
%
% Other m-files required: none
% Subfunctions: capitaliseFirstLetters, dataTypeByteSize, SDF_DataType,
%   SDF_uintX_vector, SDF_short_vector, SDF_long_vector, SDF_float_vector, SDF_double_vector,
%   SDF_DATA_SizesY, SDF_extract_labels, SDF_Y_scale.
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

numCh = length(data.SDF_CHANNEL_HDR);
%Calculate the number of vectors for this datatype
num_v = data.SDF_DATA_HDR(n_ns).total_rows * ...
    data.SDF_DATA_HDR(n_ns).total_cols;

ySize = dataTypeByteSize(data.SDF_DATA_HDR(n_ns).ydata_type);
numBytes = numCh * data.SDF_ADDTL.YDATA_block_sizes(n_ns); %Includes yPerPoint and isComplex
tmp = data.SDF_YDATA_HDR.zz_raw_data(y_idx+1:y_idx+numBytes);
% Place the relevant data into reshaped tmp array vector
tmp = reshape(tmp,...
    ySize, ...
    length(tmp)/int32(ySize));

eval(['tmp = SDF_' data.SDF_DATA_HDR(n_ns).zz_ydata_type '_vector(tmp);']);
tmp = reshape(tmp,length(tmp)/numCh, numCh);
for n_v = 1: num_v
    % Abbreviate some variable names for neatness
    dataHdr = data.SDF_DATA_HDR(n_ns);
    vctHdr = data.SDF_VECTOR_HDR(dataHdr.first_VECTOR_recordNum+n_v);
    chHdr = data.SDF_CHANNEL_HDR;
    % Define the variable name
    var_name = ['data.ext.s0d' num2str(n_ns) 'v' num2str(n_v) '.y'];
    % Scale the results
    eval([var_name ' = SDF_Y_scale(chHdr,dataHdr,vctHdr,tmp(:,n_v));']); % scale the results
    % Extract user-readable data assiciated with the processed channels/axes.
    data = SDF_extract_labels(data,var_name);
end
% Increment index
y_idx = y_idx+numBytes;
% Return values
retval = data;
end


function  [retval] = SDF_DATA_SizesY(data)
%SDF_DATA_SIZEY -Get all the Y variable and element sizes for the DATA_HRD referred data-vectors.
%  This code is an implementation of the SDF standard as defined in:
%  'User's Guide - Standard Data Format Utilities, Version B.02.01'
%  by Agilent Technologies
%  Manufacturing Part Number: 5963-1715
%  Printed in USA
%  December 1994
%  © Copyright 1989, 1991-94, 2005 Agilent Technologies, Inc.
%
% Syntax: [retval] = SDF_Multi_ChannelY(data)
%
% Inputs:
%   data - Extracted SDF headers, containing the raw y-data.
%
% Outputs:
%   retval - Extracted SDF headers, with the ollowing variables updated:
%               data.SDF_ADDTL.YDATA_var_size - Data size for the relevant data header.
%               data.SDF_ADDTL.YDATA_block_sizes - Size (in bytes) of the Y-axis Data record’s data block
%               data.SDF_ADDTL.is_scan_data - Array of 0/1 representing which DATA_HDR are scan-data
%
% Other m-files required: none
% Subfunctions: dataTypeByteSize
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

%Iterate over all the YDATA elements specified in the header
for n = 1:length(data.SDF_DATA_HDR)
    % For the present data header, calculate the data size.
    data.SDF_ADDTL.YDATA_var_size(n) = dataTypeByteSize(data.SDF_DATA_HDR(n).ydata_type);
    
    % Determine the number of Vectors for this DATA_HDR
    
    % Determine the size (in bytes) of the Y-axis Data record’s data block with
    % the following formula:
    data.SDF_ADDTL.YDATA_block_sizes(n) = ...
        int32(data.SDF_DATA_HDR(n).num_of_points) * ...
        int32(data.SDF_DATA_HDR(n).yPerPoint) * ...
        int32(data.SDF_ADDTL.YDATA_var_size(n)) * ...
        (2^(int32(data.SDF_DATA_HDR(n).yIsComplex)));
    
    % Ascertain if the data type is to be included in the scanData
    data.SDF_ADDTL.is_scan_data(n) = data.SDF_DATA_HDR(n).scanData;
end
% Return values
retval = data;
end


function  [retval] = SDF_DATA_SizesX(data)
%SDF_DATA_SIZEX - Get all the X variable and element sizes for the DATA_HRD referred data-vectors.
%  This code is an implementation of the SDF standard as defined in:
%  'User's Guide - Standard Data Format Utilities, Version B.02.01'
%  by Agilent Technologies
%  Manufacturing Part Number: 5963-1715
%  Printed in USA
%  December 1994
%  © Copyright 1989, 1991-94, 2005 Agilent Technologies, Inc.
%
% Syntax: [retval] = SDF_Multi_ChannelY(data)
%
% Inputs:
%   data - Extracted SDF headers, containing the raw y-data.
%
% Outputs:
%   retval - Extracted SDF headers, with the ollowing variables updated:
%               data.SDF_ADDTL.XDATA_var_size - Data size for the relevant data header.
%               data.SDF_ADDTL.XDATA_block_sizes - Size (in bytes) of the X-axis Data record’s data block
%               data.SDF_ADDTL.is_scan_data - Array of 0/1 representing which DATA_HDR are scan-data
%
% Other m-files required: none
% Subfunctions: dataTypeByteSize
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

% Get all the X variable and element sizes for the DATA_HRD referred data-vectors.

%Iterate over all the YDATA elements specified in the header
for n = 1:length(data.SDF_DATA_HDR)
    % For the present data header, calculate the data size.
    data.SDF_ADDTL.XDATA_var_size(n) = dataTypeByteSize(data.SDF_DATA_HDR(n).xdata_type);
    % Determine the number of Vectors for this DATA_HDR
    
    % Determine the size (in bytes) of the X-axis Data record’s data block with
    % the following formula:
    data.SDF_ADDTL.XDATA_block_sizes(n) = ...
        int32(data.SDF_DATA_HDR(n).num_of_points) * ...
        int32(data.SDF_DATA_HDR(n).xPerPoint) * ...
        int32(data.SDF_ADDTL.YDATA_var_size(n));
    
    % Ascertain if the data type is to be included in the scanData
    data.SDF_ADDTL.is_scan_data(n) = data.SDF_DATA_HDR(n).scanData;
end
% Return values
retval = data;
end

function [retval] = SDF_XDATA_Process(data,file_raw)
%SDF_XDATA_PROCESS - Identifies how x-data is stored in file_raw, and then extracts it to a
%  standardised structure.
%
%  This code is an implementation of the SDF standard as defined in:
%  'User's Guide - Standard Data Format Utilities, Version B.02.01'
%  by Agilent Technologies
%  Manufacturing Part Number: 5963-1715
%  Printed in USA
%  December 1994
%  © Copyright 1989, 1991-94, 2005 Agilent Technologies, Inc.
%
% Syntax: [retval] = SDF_XDATA_Process(data,file_raw)
%
% Inputs:
%   data - Structure containing the populated SDF header information.
%   file_raw - Raw byte array of the SDF file which is being intepreted.
%
% Outputs:
%   retval - Input variable 'data' with the x-data added
%
% Other m-files required: none
% Subfunctions: capitaliseFirstLetters, dataTypeByteSize, SDF_DataType,
%   SDF_uintX, SDF_int8, SDF_char, SDF_short, SDF_long, SDF_float, SDF_double, SDF_struct,
%   SDF_uintX_vector, SDF_short_vector, SDF_long_vector, SDF_float_vector, SDF_double_vector,
%   SDF_DATA_SizesX, SDF_extract_labels, SDF_Interpret_Template, SDF_Multi_ChannelX,
%   SDF_Multi_TraceX, SDF_Single_TraceX, SDF_X_scale.
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

% Generate array of byte-offsets for each dataset, base on DATA_HDR
data_offsets = [0 data.SDF_ADDTL.DATA_HDR_sizes(1:end-1)];

% Import the XDATA header and its raw_data based on the DATA_HDR offsets.
if data.SDF_FILE_HDR.offset_of_XDATA_record  >= 0
    data.SDF_XDATA_HDR = ...
        SDF_Interpret_Template ('SDF_XDATA_HDR', ...
        file_raw,data.SDF_FILE_HDR.offset_of_XDATA_record+1 );%+ ...
    %data_offsets(n));
    % Find the of each vector/data block so we can work at extracting the correct vectors.
    data = SDF_DATA_SizesX(data);
end

%% Include a switch statement to address the following three cases:
%      - 1 Trace (Original example), all YDATA is a single vectorize (S0D1V1)
%      - 1 Scan (No scan structure) N vectors (S0D#V#)
%      - N Scans M Vectors (S#D#V#)

% Based on the scan structure, interpred data differently.
if data.SDF_FILE_HDR.num_of_SCAN_STRUCT_record % If there is a scan structure check what type
    switch data.SDF_SCAN_STRUCT.zz_scan_type
        case 'scan'
            x_idx = 0; % Initialise the index counter for the block of data
            % Scan #, Data #, Vector #
            for n_s = 1 : data.SDF_SCAN_STRUCT.num_of_scan
                %Iterate over all the scans
                for n_d = 1 : sum(data.SDF_ADDTL.is_scan_data) %Assume sequential
                    [data,x_idx] = SDF_Multi_TraceX(data,n_s,n_d,x_idx);
                end %n_d
            end %n_s
            % Assume the non-scan data is after the scan structured data, import
            % the non-scanned data starting with the first non-scan datatype
            % Note we assume the scans are all before the non-scan.
            for n_ns = sum(data.SDF_ADDTL.is_scan_data)+1:length(data.SDF_ADDTL.is_scan_data)
                [data,x_idx] = SDF_Multi_ChannelX(data,n_ns,x_idx);
            end %n_ns
        case 'depth'
            x_idx = 0; % Initialise the index counter for the block of data
            % Data #, Scan #, Vector #
            warning('SDF_XDATA_Process: scan_type ''depth'' has not been tested.');
            for n_d = 1 : sum(data.SDF_ADDTL.is_scan_data) %Assume sequential
                %Iterate over all the scans
                for n_s = 1 : data.SDF_SCAN_STRUCT.num_of_scan
                    [data,x_idx] = SDF_Multi_TraceX(data,n_s,n_d,x_idx);
                end %n_d
            end %n_s
            % Assume the non-scan data is after the scan structured data, import
            % the non-scanned data starting with the first non-scan datatype
            % Note we assume the scans are all before the non-scan.
            for n_ns = sum(data.SDF_ADDTL.is_scan_data)+1:length(data.SDF_ADDTL.is_scan_data)
                [data,x_idx] = SDF_Multi_ChannelX(data,n_ns,x_idx);
            end %n_ns
        otherwise
            % Do nothing
    end % switch data.SDF_SCAN_STRUCT.zz_scan_type
else % There is no scan struct, read y trace
    data = SDF_Single_TraceX(data); % Old way for single trace
end % if data.SDF_FILE_HDR.num_of_SCAN_STRUCT_record
retval = data;
end



function [retval,x_idx] = SDF_Multi_TraceX(data,n_s,n_d,x_idx)
%SDF_MULTI_TRACEX - Reconstructs the scan X-data for the datasets which are part of a multiscan file.
%  This code is an implementation of the SDF standard as defined in:
%  'User's Guide - Standard Data Format Utilities, Version B.02.01'
%  by Agilent Technologies
%  Manufacturing Part Number: 5963-1715
%  Printed in USA
%  December 1994
%  © Copyright 1989, 1991-94, 2005 Agilent Technologies, Inc.
%
% Syntax: [retval] = SDF_Multi_TraceX(data)
%
% Inputs:
%   data - Extracted SDF headers, containing the raw x-data.
%   n_s - Scan number (must be greater than zero)
%   n_d - Data header number (must be greater than zero)
%   y_idx - index of the first byte relevant to the Multi Trace input data.SDF_YDATA_HDR.zz_raw_data
%
% Outputs:
%   retval - Extracted SDF headers, and processed x-axis traces.
%
% Other m-files required: none
% Subfunctions: capitaliseFirstLetters, dataTypeByteSize, SDF_DataType, SDF_uintX_vector,
% SDF_short_vector, SDF_long_vector, SDF_float_vector, SDF_double_vector, SDF_extract_labels.
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

%Calculate the number of vectors for this datatype
num_v = data.SDF_DATA_HDR(n_d).total_rows * ...
    data.SDF_DATA_HDR(n_d).total_cols;

switch data.SDF_DATA_HDR(n_d).zz_xResolution_type
    case {'linear'}
        n = 0:1.0:data.SDF_DATA_HDR(n_d).num_of_points-1;
        data.SDF_XDATA_HDR(n_d).zz_data = ...
            data.SDF_DATA_HDR(n_d).abscissa_firstX + ...
            data.SDF_DATA_HDR(n_d).abscissa_deltaX * n;
        for n_v = 1:num_v
            % Define the variable name
            var_name = ['data.ext.s' num2str(n_s) ...
                'd' num2str(n_d) 'v' num2str(n_v) '.x'];
            eval([var_name ...
                ' = data.SDF_XDATA_HDR(n_d).zz_data;']); % import the results
            % Extract user-readable data assiciated with the processed channels/axes.
            data = SDF_extract_labels(data,var_name);
        end
    case {'logarithmic'}
        n = 0:data.SDF_XDATA_HDR.num_of_points-1;
        data.SDF_XDATA_HDR(n_d).zz_data = ...
            data.SDF_DATA_HDR(n_d).abscissa_firstX * ...
            data.SDF_DATA_HDR(n_d).abscissa_deltaX ^ n;
        for n_v = 1:num_v
            % Define the variable name
            var_name = ['data.ext.s' num2str(n_s) ...
                'd' num2str(n_d) 'v' num2str(n_v) '.x'];
            eval([var_name ...
                ' = data.SDF_XDATA_HDR(n_d).zz_data;']); % import the results
            % Extract user-readable data assiciated with the processed channels/axes.
            data = SDF_extract_labels(data,var_name);
        end
    otherwise
        % Read the block size we have calculated
        blk_size = data.SDF_ADDTL.XDATA_block_sizes(n_d);
        % Import all the vectors reshaped to be in rows
        for n_v = 1:num_v
            strt = x_idx+(n_v-1) * blk_size+1;
            stp = x_idx+(n_v) * blk_size;
            % Define the variable name
            var_name = ['data.ext.s' num2str(n_s) ...
                'd' num2str(n_d) 'v' num2str(n_v) '.x'];
            dataHdr = data.SDF_DATA_HDR(n_d);
            vctHdr = data.SDF_VECTOR_HDR(dataHdr.first_VECTOR_recordNum+n_v);
            chHdr(1) = data.SDF_CHANNEL_HDR(vctHdr.the_CHANNEL_record(1)+1);
            if (vctHdr.the_CHANNEL_record(2) > -1)
                chHdr(2) = data.SDF_CHANNEL_HDR(vctHdr.the_CHANNEL_record(2)+1);
            end
            % Place the relevant data into reshaped tmp array vector
            tmp = reshape(data.SDF_XDATA_HDR.zz_raw_data(strt:stp),...
                data.SDF_ADDTL.XDATA_var_size(n_d), ...
                blk_size/int32(data.SDF_ADDTL.XDATA_var_size(n_d)));
            % Convert the vector to datatype and scale in single command.
            eval([var_name ' = SDF_' data.SDF_DATA_HDR(n_d).zz_xdata_type '_vector(tmp);']); % import the results
            % Extract user-readable data assiciated with the processed channels/axes.
            data = SDF_extract_labels(data,var_name);
        end
        x_idx = stp;
end
% Increment the index counter
retval = data;
end

function [retval,x_idx] = SDF_Multi_ChannelX(data,n_ns,x_idx)
%SDF_MULTI_CHANNELX - Reconstructs non-scan X-data for the datasets which are part of a multiscan file.
%  This code is an implementation of the SDF standard as defined in:
%  'User's Guide - Standard Data Format Utilities, Version B.02.01'
%  by Agilent Technologies
%  Manufacturing Part Number: 5963-1715
%  Printed in USA
%  December 1994
%  © Copyright 1989, 1991-94, 2005 Agilent Technologies, Inc.
%
% Syntax: [retval] = SDF_Multi_ChannelX(data)
%
% Inputs:
%   data - Extracted SDF headers, containing the raw y-data.
%   n_ns - Data header number (must be greater than zero)
%   x_idx - index of the first byte relevant to the Multi Channel input data.SDF_XDATA_HDR.zz_raw_data
%
% Outputs:
%   retval - Extracted SDF headers, and processed y-axis traces.
%
% Other m-files required: none
% Subfunctions: capitaliseFirstLetters, dataTypeByteSize, SDF_DataType, SDF_uintX_vector,
%   SDF_short_vector, SDF_long_vector, SDF_float_vector, SDF_double_vector, SDF_DATA_SizesX,
%   SDF_extract_labels, SDF_X_scale.
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

%Calculate the number of vectors for this datatype
num_v = data.SDF_DATA_HDR(n_ns).total_rows * data.SDF_DATA_HDR(n_ns).total_cols;

switch data.SDF_DATA_HDR(n_ns).zz_xResolution_type
    case {'linear'}
        n = 0:data.SDF_DATA_HDR(n_ns).num_of_points-1;
        data.SDF_XDATA_HDR(n_ns).zz_data = ...
            data.SDF_DATA_HDR(n_ns).abscissa_firstX + ...
            data.SDF_DATA_HDR(n_ns).abscissa_deltaX * n;
        for n_v = 1:num_v
            % Define the variable name
            var_name = ['data.ext.s0' ...
                'd' num2str(n_ns) 'v' num2str(n_v) '.x'];
            eval([var_name ...
                ' = data.SDF_XDATA_HDR(n_ns).zz_data;']); % import the results
            % Extract user-readable data assiciated with the processed channels/axes.
            data = SDF_extract_labels(data,var_name);
        end
    case {'logarithmic'}
        n = 0:data.SDF_XDATA_HDR.num_of_points-1;
        data.SDF_XDATA_HDR(n_ns).zz_data = ...
            data.SDF_DATA_HDR(n_ns).abscissa_firstX * ...
            data.SDF_DATA_HDR(n_ns).abscissa_deltaX ^ n;
        for n_v = 1:num_v
            % Define the variable name
            var_name = ['data.ext.s0' ...
                'd' num2str(n_ns) 'v' num2str(n_v) '.x'];
            eval([var_name ...
                ' = data.SDF_XDATA_HDR(n_ns).zz_data;']); % import the results
            % Extract user-readable data assiciated with the processed channels/axes.
            data = SDF_extract_labels(data,var_name);
        end
    otherwise
        % Issue warning since this section of code was never tested (lack of test input files)
        warning([data.SDF_ADDTL.filename ': An untested section of code was executed in '...
            'SDF_Multi_ChannelX: data.SDF_DATA_HDR(' num2str(n_ns) ').zz_xResolution_type = ''' ...
            data.SDF_DATA_HDR(n_ns).zz_xResolution_type '''']);
        % Account for the byte Size of the variable
        xSize = dataTypeByteSize(data.SDF_DATA_HDR(n_ns).xdata_type);
        
        % Calculate the number of bytes to read
        numBytes = int32(data.SDF_DATA_HDR(n_ns).num_of_points * xSize);
        % Get the data
        tmp = data.SDF_YDATA_HDR.zz_raw_data(x_idx+1:x_idx+numBytes);
        % Place the relevant data into reshaped tmp array vector
        tmp = reshape(tmp,xSize,length(tmp)/int32(xSize));
        for n_v = 1:num_v
            % Define the variable name path to the standardise export location
            var_name = ['data.ext.s0' 'd' num2str(n_ns) 'v' num2str(n_v) '.x'];
            %Convert the data to the correct data type
            eval([var_name ' = SDF_' data.SDF_DATA_HDR(n_ns).zz_xdata_type '_vector(tmp);']);
            % Extract user-readable data assiciated with the processed channels/axes.
            data = SDF_extract_labels(data,var_name);
            % Increment to the next block of x_idx
            x_idx = x_idx+numBytes;
        end
end
retval = data; % Return data
end

function [retval] = SDF_Single_TraceX(data)
%SDF_SINGLE_TRACEX - Reconstructs the X-data for the simplest type of SDF file, which contains a
%  single trace.
%  This code is an implementation of the SDF standard as defined in:
%  'User's Guide - Standard Data Format Utilities, Version B.02.01'
%  by Agilent Technologies
%  Manufacturing Part Number: 5963-1715
%  Printed in USA
%  December 1994
%  © Copyright 1989, 1991-94, 2005 Agilent Technologies, Inc.
%
% Syntax: [retval] = SDF_Single_TraceX(data)
%
% Inputs:
%   data - Extracted SDF headers, containing the raw x-data.
%   file_raw - Raw byte array of the SDF file which is being intepreted.
%
% Outputs:
%   retval - Extracted SDF headers, and processed x-axis traces.
%
% Other m-files required: none
% Subfunctions: capitaliseFirstLetters, dataTypeByteSize, SDF_DataType,
%   SDF_uintX_vector, SDF_short_vector, SDF_long_vector, SDF_float_vector, SDF_double_vector,
%   SDF_DATA_SizesX, SDF_extract_labels.
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence


% Define the variable name path to the standardise export location
var_name = 'data.ext.s0d1v1.x';

switch data.SDF_DATA_HDR.zz_xResolution_type
    case {'linear'} % X-axis has a linear scale
        n = 0:data.SDF_DATA_HDR.num_of_points-1;
        data.SDF_XDATA_HDR.zz_data = ...
            data.SDF_DATA_HDR.abscissa_firstX + ...
            data.SDF_DATA_HDR.abscissa_deltaX * n;
    case {'logarithmic'} % X-axis has a log scale
        n = 0:data.SDF_XDATA_HDR.num_of_points-1;
        data.SDF_XDATA_HDR.zz_data = ...
            data.SDF_DATA_HDR.abscissa_firstX * ...
            data.SDF_DATA_HDR.abscissa_deltaX ^ n;
        %    case {'arbitrary, one per file'}
        %    case {'arbitrary, one per data type'}
        %    case {'arbitrary, one per trace'}
    otherwise
        % Issue warning since this section of code was never tested (lack of test input files)
        warning([data.SDF_ADDTL.filename ': An untested section of code was executed in '...
            'SDF_Single_TraceX: data.SDF_DATA_HDR.zz_xResolution_type = ''' ...
            data.SDF_DATA_HDR.zz_xResolution_type '''']);
        % Issue warning since this section is not fully implemented
        warning([data.SDF_ADDTL.filename ': An un-implmented section of code was executed in '...
            'SDF_Single_TraceX: data.SDF_DATA_HDR.zz_xResolution_type = ''' ...
            data.SDF_DATA_HDR.zz_xResolution_type '''']);
        
        % Calculate how many bytes there are per datapoint
        var_size = dataTypeByteSize(data.SDF_DATA_HDR.xdata_type);
        
        % Reshape the _raw_data to make a 2D array of bytes relevant to x-axis
        data.SDF_XDATA_HDR.zz_raw_data = ...
            reshape(data.SDF_XDATA_HDR.zz_raw_data, ...
            var_size,length(data.SDF_XDATA_HDR.zz_raw_data)/var_size);
        
        % Convert the smaller arrays to their appropriate data types.
        for n = 1:size(data.SDF_XDATA_HDR.zz_raw_data,2)
            eval(['data.SDF_XDATA_HDR.zz_data(n) = ' ...
                'SDF_' data.SDF_DATA_HDR.zz_xdata_type ...
                '(data.SDF_XDATA_HDR.zz_raw_data(:,n));']);
        end
        
        % Reshape the array to de-interlieve the datasets.
        data.SDF_XDATA_HDR.zz_data = ...
            reshape(data.SDF_XDATA_HDR.zz_data, ...
            data.SDF_DATA_HDR.xPerPoint, ...
            length(data.SDF_XDATA_HDR.zz_data)/data.SDF_DATA_HDR.xPerPoint);
end
% Save the data to the standardise export location
eval([var_name ' = data.SDF_XDATA_HDR.zz_data;']);
% Extract user-readable data assiciated with the processed channels/axes.
data = SDF_extract_labels(data,var_name);
retval = data; % Return data
end

function [retval] = SDF_extract_labels(data,var_name)
%SDF_EXTRACT_LABELS - Extracts the trace units (e.g. V^2/Hz), and other relevant descriptors, such
% as channel/module identifiers for the trace specified in 'var_name'.
%  This code is an implementation of the SDF standard as defined in:
%  'User's Guide - Standard Data Format Utilities, Version B.02.01'
%  by Agilent Technologies
%  Manufacturing Part Number: 5963-1715
%  Printed in USA
%  December 1994
%  © Copyright 1989, 1991-94, 2005 Agilent Technologies, Inc.
%
% Syntax: [retval] = SDF_extract_labels(data,var_name)
%
% Inputs:
%   data - Extracted SDF headers, containing the headers, and extracted trace.
%   var_name - Trace structure variable name, matching the format 'data.ext.s#d#v#.axis', where for:
%              s (scan number), # >= 0, since 'Scan 0' represents non-scan data, e.g. 's0'
%              d (data header number), # > 0, since it is an index, e.g. 'd2'
%              v (vector header number), # > 0, since it is an index, e.g. 'v2'
%              axis - specifies which axis is beinf processed: 'x' or 'y'
%
% Outputs:
%   retval - Extracted SDF headers, and processed labels.
%            The extracted data is stores in retval.ext.s#d#v# as follows:
%              - xUnit, the unit-label associated with the x-axis by the SDF file
%              - xResolution_type, identifies if data is linear, logarithmic or one of
%                  the three arbitrary formats.
%              - xLabel, the domain of the trace: 'Time', 'Frequency' or blank.
%                  (Based on the x-channel units)
%              - ch#num, the physical channel number from which the trace was taken (# is 1 or 2)
%              - ch#label, the stored label of the physical channel from which the trace was taken.
%              - ch#moduleId, location of the channel in the instrument (or instrument i.d.).
%              - ch#serialNum, instrument or module serial number.
%              - yUnit, the unit-label associated with the y-axis by the SDF file.
%              - yLabel, DATA_HDR dataType variable with the first letter of each word capitalised.
%              - dataType, the DATA_HDR dataType variable.
%              - dataTitle, the DATA_HDR dataTitle variable.
%              - yIsPowerData, the DATA_HDR yIsPowerData variable.
%              - idx, the vector of valid data identified by the instrument:
%                  [SDF_MEAS_HDR.startFreqIndex:SDF_MEAS_HDR.stopFreqIndex)+1]
%              - DATA_HDRnum, the DATA_HDR index (for convenience)
%              - VECTOR_HDRnum, the VECTOR_HDR index (for convenience)
%              - SCANnum, the SCAN number (for convenience)
%              - traceList, a list of every trace and its associated dataType.
%
% Other m-files required: none
% Subfunctions: capitaliseFirstLetters.
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

% Interpret the variable name to extract the scan, data and vector numbers.
[s] = sscanf(var_name(1:end-2),'data.ext.s%dd%dv%d');
d = s(2); %DATA_HDR #
v = s(3); %VECTOR_HDR #
s = s(1); %SCAN #
x_or_y = var_name(end); %Determine if we are labelling the x or y axis

% Allow function to operate differently if x or y dataset is being processed
% This accounts for the fact that there is no x data extracted until the endof the program
switch x_or_y
    case 'x'
        % xLabel is always defined
        eval([var_name(1:end-1) 'xUnit = data.SDF_DATA_HDR(d).xUnit.label;']);
        % Save the resolution type (linear or logarithmic)
        eval([var_name(1:end-1) ...
            'xResolution_type = data.SDF_DATA_HDR(d).zz_xResolution_type;']);
        % Use the x label to identify is we are in the frequency or time domain
        if strcmp(data.SDF_DATA_HDR(d).xUnit.label,'Hz') % Frequency domain
            eval([var_name(1:end-1) 'xLabel = ''Frequency'';']);
        elseif strcmp(data.SDF_DATA_HDR(d).xUnit.label,'s') % Time domain
            eval([var_name(1:end-1) 'xLabel = ''Time'';']);
        else % Unknown unit used
            % Provide a warning
            warning([data.SDF_ADDTL.filename ': An unknown x-axis unit ''' ...
                data.SDF_DATA_HDR(d).xUnit.label ''' was used in ' ...
                'data.SDF_DATA_HDR(' num2str(d) ').xUnit.label, ' var_name(1:end-1) 'xLabel' ...
                ' has been left blank.']);
            % Provide blank variable
            eval([var_name(1:end-1) 'xLabel = '''';']);
        end
    case 'y'
        % Save the instrument dependent channel Label, model and serial number
        flag = [0 0]; % A flag to identify which channels are used
        for ch = 1 : 2 % Iterate over numerator and denominator channels
            tmp = data.SDF_VECTOR_HDR(v).the_CHANNEL_record(ch)+1; %Included 0 offset
            if tmp > 0
                % Channel Number
                eval([var_name(1:end-1) 'ch' num2str(ch) 'num = tmp;']);
                % Channel label
                eval([var_name(1:end-1) 'ch' num2str(ch) ...
                    'label = data.SDF_CHANNEL_HDR(tmp).channelLabel;']);
                % Module ID
                eval([var_name(1:end-1) 'ch' num2str(ch) ...
                    'moduleId = data.SDF_CHANNEL_HDR(tmp).moduleId;']);
                % Serial number
                eval([var_name(1:end-1) 'ch' num2str(ch) ...
                    'serialNum = data.SDF_CHANNEL_HDR(tmp).serialNum;']);
            else % Flag is the channel is unused
                flag(ch) = 1;
            end
        end
        
        % Check for presense of a yUnit label
        if length(data.SDF_DATA_HDR(d).yUnit.label)<1
            flag=[2 1]*flag';
            switch flag % If no units have been upplied, make them up using the channel data.
                case 0%[0 0] % Both channels are used (E.g. A transfer function)
                    eval([var_name(1:end-1) 'yUnit = [data.SDF_CHANNEL_HDR(1).engUnit.label ''/''  data.SDF_CHANNEL_HDR(2).engUnit.label];']);
                case 1%[0 1] % Only denominator channel is used (unlikely to ever occur)
                    eval([var_name(1:end-1) 'yUnit = data.SDF_CHANNEL_HDR(2).engUnit.label;']);
                case 2%[1 0] % Only numerator channel is used (E.g. noise/time measurement)
                    eval([var_name(1:end-1) 'yUnit = data.SDF_CHANNEL_HDR(1).engUnit.label;']);
                case 3%[1 1] % Neither channel is used (unlikely to ever occur, used as a catchall)
                    eval([var_name(1:end-1) 'yUnit = data.SDF_DATA_HDR(d).yUnit.label;']);
            end
        else % Use the yUnit label if it has been supplied
            eval([var_name(1:end-1) 'yUnit = data.SDF_DATA_HDR(d).yUnit.label;']);
        end
        % Use the datatype as a temporary label for the dataset.
        eval([var_name(1:end-1) 'yLabel = capitaliseFirstLetters(data.SDF_DATA_HDR(d).zz_dataType);']);
        % Copy the datatype without capitalisation
        eval([var_name(1:end-1) 'dataType = data.SDF_DATA_HDR(d).zz_dataType;']);
        % Copy the userdefined title for the dataset
        eval([var_name(1:end-1) 'dataTitle = data.SDF_DATA_HDR(d).dataTitle;']);
        % Record if the data is power data
        eval([var_name(1:end-1) 'yIsPowerData = data.SDF_DATA_HDR(d).zz_yIsPowerData;']);
        % Record the valid indexes
        eval([var_name(1:end-1) ...
            'idx=(data.SDF_MEAS_HDR.startFreqIndex:data.SDF_MEAS_HDR.stopFreqIndex)+1;']);
        % Save the DATA_HDR index (for convenience)
        eval([var_name(1:end-1) 'DATA_HDRnum = d;']);
        % Save the VECTOR_HDR index (for convenience)
        eval([var_name(1:end-1) 'VECTOR_HDRnum = v;']);
        % Save the SCAN number (for convenience)
        eval([var_name(1:end-1) 'SCANnum = s;']);
        % Save the tracetype and varname to a list of cells in data.ext
        data.ext.traceList = [data.ext.traceList ;...
            {['s',num2str(s, '%d'),'d',num2str(d, '%d'),'v',num2str(v,'%d')]},...
            { data.SDF_DATA_HDR(d).zz_dataType}];
end
retval = data;
end

function txt = capitaliseFirstLetters(txt)
%CAPITALISEFIRSTLETTERS - This function capitalises the first character in every word of a text string.
%
% Syntax: [txt] = capitaliseFirstLetters(txt)
%
% Inputs:
%   txt - Test requiring reformatting
%
% Outputs:
%   txt - The processed text string.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

txt(1) = upper(txt(1)); % Always capitalise the first char in the string_fill_char
idx = strfind(txt,' ')+1; % Identify first letter after each space
idx = idx(idx<length(txt)); % Remove case where last character is a space
for n = idx % Iterate over the characters to capitalise
    txt(n) = upper(txt(n));
end
end

function varSize = dataTypeByteSize(dataType)
%DATATYPEBYTESIZE - Calculates the number of bytes associated with dataType.
% dataType = 1, short  (varSize = 2 bytes)
% dataType = 2, long   (varSize = 4 bytes)
% dataType = 3, float  (varSize = 4 bytes)
% dataType = 4, double (varSize = 8 bytes)
%
% Syntax: [varSize] = dataTypeByteSize(dataType)
%
% Inputs:
%   dataType - Integer of value 1 to 4 (see SDF standard)
% Outputs:
%   varSize - Positive integer representing the number of bytes required.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 04/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

if dataType == 3 % 3 = float
    varSize = 4;
elseif(dataType  <= 4) % 1 = short, 2 = long, 4 = double
    varSize = dataType * 2;
else % Other datatypes have not been defined in the SDF standard.
    error('Unknown dataType');
end
end

function [retval] = SDF_Y_scale(chHdr,dataHdr,vctHdr,y)
%SDF_Y_SCALE - Scales and re-arranges trace to the appropriate units/shape.
%
% Syntax: [retval] = SDF_Y_scale(chHdr,dataHdr,vctHdr,y)
%
% Inputs:
%   chHdr - Channel headers relevant to the data.
%   dataHdr - Data header relevant to the data.
%   vctHdr - Vector header relevant to the data.
%   y - Unscaled data
% Outputs:
%   retval - Scaled and re-arranged results.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: justin.dinale@dsto.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 06/04/2017
% © Copyright 2016-2017 Commonwealth of Australia, represented by the Department of Defence

%SDF_Y_scale
% Calculate scaling factor base on dataheader type
y = double(y);

switch dataHdr.zz_dataType
    case 'time'
        % Since only one channel is used, identify which is valid
        ch_idx = vctHdr.the_CHANNEL_record > -1;
        % Extract relevant channel offset and scale, apply them
        retval = chHdr(ch_idx).channelOffset + y * chHdr(ch_idx).channelScale;
    case 'overload data'
        retval = y;
    case 'decimated time data'
        retval.overloadFlag = y(5:5:end);
        % the minimume/maximum values are in integer for, so need scaling.
        % Since only one channel is used, identify which is valid
        ch_idx = vctHdr.the_CHANNEL_record > -1;
        % Extract relevant channel offset and scale, apply them
        y = chHdr(ch_idx).channelOffset + y * chHdr(ch_idx).channelScale;
        retval.decimatedRealValue = y(1:5:end);
        retval.decimatedImaginaryValue = y(2:5:end);
        warning(['Decimated time data has been extracted. '...
            'This code was not tested. Please verify.']);
    case 'compressed time data'
        retval.overloadFlag = y(5:5:end);
        % the minimume/maximum values are in integer for, so need scaling.
        % Since only one channel is used, identify which is valid
        ch_idx = vctHdr.the_CHANNEL_record > -1;
        % Extract relevant channel offset and scale, apply them
        y = chHdr(ch_idx).channelOffset + y * chHdr(ch_idx).channelScale;
        retval.minimumRealValue = y(1:5:end);
        retval.minimumImaginaryValue = y(2:5:end);
        retval.maximumRealValue = y(3:5:end);
        retval.maximumImaginaryValue = y(4:5:end);
    case 'tachometer data'
        retval.number_of_tach_points = (dataHdr.total_rows * dataHdr.total_cols) * ...
            dataHdr.num_of_points + dataHdr.last_valid_index;
        retval.tach_pulses_per_rev = dataHdr.abscissa_deltaX;
        retval.tach_frequency = dataHdr.abscissa_firstX;
        retval.tach_time_sec = y/retval.tach_frequency;
        retval.tach_pulse_delta_sec = diff(retval.tach_time_sec);
        retval.tach_pulses_per_sec = 1./retval.tach_pulse_delta_sec;
        retval.tach_pulses_per_min = 60 * retval.tach_pulses_per_sec;
        retval.tach_rpm = retval.tach_pulses_per_min / retval.tach_pulses_per_rev;
        retval.userDelay = chHdr.userDelay; % Delay between tachometer zero count and start of capture.
        warning(['Tachometer data has been extracted. ' ...
            'This code was not tested. Please verify.']);
    case 'external trigger data'
        retval.number_of_ext_trigger_points = (dataHdr.total_rows * dataHdr.total_cols) * ...
            dataHdr.num_of_points + dataHdr.last_valid_index;
        retval.ext_trigger_pulses_per_rev = dataHdr.abscissa_deltaX;
        retval.ext_trigger_frequency = dataHdr.abscissa_firstX;
        retval.ext_trigger_time_sec = y/retval.ext_trigger_frequency;
        retval.ext_trigger_pulse_delta_sec = diff(retval.ext_trigger_time_sec);
        retval.ext_trigger_pulses_per_sec = 1./retval.ext_trigger_pulse_delta_sec;
        retval.ext_trigger_pulses_per_min = 60 * retval.ext_trigger_pulses_per_sec;
        retval.ext_trigger_rpm = retval.ext_trigger_pulses_per_min / retval.ext_trigger_pulses_per_rev;
        retval.userDelay = chHdr.userDelay; % Delay between tachometer zero count and start of capture.
        warning(['External trigger data has been extracted. ' ...
            'This code was not tested. Please verify.']);
        
    otherwise % Assume everything else is treated the same way as frequency domain
        %% The following window correction will require further work. The documentation
        % is unclear  on how to identify which type of correction is required when a
        % correction has not already been applied.
        % Correction factors
        WINDOW_corr(1:2) = 1; % Initialise flag variable
        if  or(strcmp(dataHdr.zz_domain,'frequency') , ...
                strcmp(dataHdr.zz_domain,'order'))
            for m = (1:2) %Iterate over numberator and denominator.
                switch chHdr(m).window.windowCorrMode
                    case 0 % No window
                        WINDOW_corr(m) = 1;
                    case 1 % Narrow band window
                        WINDOW_corr(m) = chHdr(m).window.narrowBandCorr;
                    case 2 % Broadband window
                        WINDOW_corr(m) = chHdr(m).window.wideBandCorr;
                end %switch chHdr(m).window.windowCorrMode
            end % for m = (1:2)
        end %if  or(strcmp(data.SDF_DATA_HDR(n).zz_domain,'frequency') , ...
        %    strcmp(data.SDF_DATA_HDR(n).zz_domain,'order'))
        
        % Calculare the integer to engineerig units and 'power' of the channel
        % corrections
        EU_corr(1) = chHdr(1).int2engrUnit;
        EU_corr(2) = chHdr(2).int2engrUnit;
        % Divide VECTOR_HDR.powrOfChan by 48 as specified in the standard
        POW_corr(1) = double(vctHdr.powrOfChan(1)/48);
        POW_corr(2) = double(vctHdr.powrOfChan(2)/48);
        if vctHdr.the_CHANNEL_record(2) == -1
            % Forcibly disable WINDOW correction, due to them already having been applied.
            % This hack is required when only one channel is used. See page B-30 (Note & final paragraph)
            WINDOW_corr = [1 1];
            if 0, warning(['The property ''windowCorrMode'' has been ignored in ' ...
                'SDF_Single_TraceY. (See B-30 of SDF Standard)']); end
        end
        % Make array of 1's and 0's on whether to apply correction.
        use_channel = (vctHdr.the_CHANNEL_record ~= -1);
        % Page B-33 '4. Create a correction factor by multiplying the two chanenl correction factors'
        correction = prod((WINDOW_corr./EU_corr).^(use_channel .* POW_corr));
        % Page B-31 disagrees with B-33 (unless reciproal is part of SDF_VECTOR_HDR.powrOfChan(2))
        %correction = (WINDOW_corr./EU_corr).^(use_channel .* POW_corr);
        %correction = correction(1)/correction(2);
        
        % If data is complex then take real and imaginary parts.
        if dataHdr.yIsComplex
            y = y(1:2:end-1) + 1i * y(2:2:end);
        end
        % Reshape the data to get the multiple traces as specified in
        % SDF_DATA_HDR(n).yPerPoint
        y = reshape(y, dataHdr.yPerPoint, length(y)/dataHdr.yPerPoint);
        % Apply the correction
        retval = y * correction;
        
end
end

function [retval] = setGlobal(var,var1)
% Debugging function.  Push a local variale to be global, allowing it to be accessed.
global debug_var
global debug_var1

debug_var = var;
debug_var1 = var1;
end

function [retval] = SDF_template (template)
%SDF_TEMPLATE - Provides an empty template containing the Standard Data Format (SDF) header 
%  information required to read generic HP/Agilent/Keysight SDF files.
%  The header structure information is detailed in:
%  "User's Guide - Standard Data Format Utilities, Version B.02.01"
%  by Agilent Technologies 
%  Manufacturing Part Number: 5963-1715
%  Printed in USA
%  December 1994
%  © Copyright 1989, 1991-94, 2005 Agilent Technologies, Inc.
% 
% Syntax: [retval] = SDF_template (template)
%
% Inputs:
%   template - String of the desired template, options are:
%                 - SDF_FILE_HDR:    Provides an index to the file.
%                 - SDF_MEAS_HDR:    Contains settings of measurement parameters.
%                 - SDF_DATA_HDR:    Tells you how reconstruct one block of measurement results 
%                                    (x- and y-axis values for every point of every trace).
%                 - SDF_VECTOR_HDR:  Tells you which channel (or pair of channels) provided data 
%                                    for a single trace.
%                 - SDF_CHANNEL_HDR: Contains channel-specific information for one channel 
%                                    used in the measurement.
%                 - SDF_SCAN_STRUCT: Tells you how vectors are organized in the Y-axis Data 
%                                    record when the measurement includes multiple scans of data.
%                 - SDF_SCANS_BIG:   Extended scan header which tells how vectors are organized 
%                                    in the Y-axis data record when the measurement may include 
%                                    more than 32767 scans of data.
%                 - SDF_SCANS_VAR:   Contains a number that identifies every scan 
%                                    (scan time, RPM, or scan number).
%                 - SDF_COMMENT_HDR: Contains text information associated with the file.
%                 - SDF_XDATA_HDR:   Contains the x-axis data needed to reconstruct any trace.
%                 - SDF_YDATA_HDR:   Contains the y-axis data needed to reconstruct any trace.
%                 - SDF_UNIQUE:      Makes the SDF file format flexible. The eight common record
%                                    types define parameters that are common to many instruments 
%                                    and systems. However, a particular instrument or system may 
%                                    need to save and recall the states of additional parameters. 
%                                    These states can reside in Unique records.
%                 - SDF_UNIT:        Contains eng. units & scaling information for the traces.
%                 - SDF_WINDOW:      Contains windowing information for frequency domain traces.
% Outputs:
%   data - Extracted SDF headers, and processed signal traces.
%   file_raw - Unprocessed binary stream.
%
% Other m-files required: none
% Subfunctions: SDF_CHANNEL_HDR_Template, SDF_COMMENT_HDR_Template, SDF_DATA_HDR_Template,
%   SDF_FILE_HDR_Template, SDF_MEAS_HDR_Template, SDF_SCAN_STRUCT_Template, SDF_SCANS_BIG_Template,
%   SDF_SCANS_VAR_Template, SDF_UNIT_Template, SDF_VECTOR_HDR_Template, SDF_WINDOW_Template,
%   SDF_XDATA_HDR_Template, SDF_YDATA_HDR_Template, SDF_UNIQUE_Template.
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: N/A
% © Copyright 2016 Commonwealth of Australia, represented by the Department of Defence

%   This file is part of SDF Importer.
% 
% Copyright 2016-2018 Commonwealth of Australia, represented by the Department of Defence.
%
% Redistribution and use in source and binary forms, with or without modification, are 
% permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this list of 
%    conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice, this list 
%    of conditions and the following disclaimer in the documentation and/or other materials 
%    provided with the distribution.
%
% 3. Neither the name of the copyright holder nor the names of its contributors may be used 
%    to endorse or promote products derived from this software without specific prior 
%    written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY 
% AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% Wrapper function to call any of the Template functions within the file.
  eval(['retval=' template '_Template();']); 
end

function  SDF_CHANNEL_HDR=SDF_CHANNEL_HDR_Template()
%SDF_CHANNEL_HDR_TEMPLATE - Provides a empty template containing the requested header information
%
%                 - SDF_CHANNEL_HDR: Contains channel-specific information for one channel 
%                                    used in the measurement.
% 
% Syntax: [retval] = SDF_CHANNEL_HDR_TEMPLATE()
%
% Inputs: none
%
% Outputs:
%   retval - Extracted SDF header information.
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: 16/03/2018
% © Copyright 2016-2018 Commonwealth of Australia, represented by the Department of Defence
%
% v1.0   27/04/2018 Changed version number to match MATLAB Central
% v0.2.3 16/03/2018 Changed SDF_CHANNEL_HDR.overloaded template entry to
% deal with error in an overloaded file. The value '3' (out of 0-1 range)
% was presen for a file which exhibited 'overloaded' on the instrument
% screen.

SDF_CHANNEL_HDR(1).FieldIndex=1;
SDF_CHANNEL_HDR(1).BinaryIndex=[1 2];
SDF_CHANNEL_HDR(1).FieldName='recordType';
SDF_CHANNEL_HDR(1).DataType='short';
SDF_CHANNEL_HDR(1).RangeUnits='14';

SDF_CHANNEL_HDR(2).FieldIndex=2;
SDF_CHANNEL_HDR(2).BinaryIndex=[3 6];
SDF_CHANNEL_HDR(2).FieldName='recordSize';
SDF_CHANNEL_HDR(2).DataType='long';
SDF_CHANNEL_HDR(2).RangeUnits='212 bytes';

SDF_CHANNEL_HDR(3).FieldIndex=3;
SDF_CHANNEL_HDR(3).BinaryIndex=[7 10];
SDF_CHANNEL_HDR(3).FieldName='unique_record';
SDF_CHANNEL_HDR(3).DataType='long';
SDF_CHANNEL_HDR(3).RangeUnits='-1:2^31-1';
SDF_CHANNEL_HDR(3).Description=...
['byte offset from the beginning of the file to ' ...
'a record containing an instrument-specific vector header. ' ...
'May be ignored if recalled into a different type instrument.'];

SDF_CHANNEL_HDR(4).FieldIndex=4;
SDF_CHANNEL_HDR(4).BinaryIndex=[11 40];
SDF_CHANNEL_HDR(4).FieldName='channelLabel';
SDF_CHANNEL_HDR(4).DataType='char[30]';
SDF_CHANNEL_HDR(4).RangeUnits='i.d.';
SDF_CHANNEL_HDR(4).Description=['channel documentation'];

SDF_CHANNEL_HDR(5).FieldIndex=5;
SDF_CHANNEL_HDR(5).BinaryIndex=[41 52];
SDF_CHANNEL_HDR(5).FieldName='moduleId';
SDF_CHANNEL_HDR(5).DataType='char[12]';
SDF_CHANNEL_HDR(5).RangeUnits='';
SDF_CHANNEL_HDR(5).Description=['location of channel in instrument.'];

SDF_CHANNEL_HDR(6).FieldIndex=6;
SDF_CHANNEL_HDR(6).BinaryIndex=[53 64];
SDF_CHANNEL_HDR(6).FieldName='serialNum';
SDF_CHANNEL_HDR(6).DataType='char[12]';
SDF_CHANNEL_HDR(6).RangeUnits='';
SDF_CHANNEL_HDR(6).Description=['instrument (or module) serial number.'];

SDF_CHANNEL_HDR(7).FieldIndex=7;
SDF_CHANNEL_HDR(7).BinaryIndex=[65 88];
SDF_CHANNEL_HDR(7).FieldName='window';
SDF_CHANNEL_HDR(7).DataType='struct';
SDF_CHANNEL_HDR(7).RangeUnits='see SDF_WINDOW';
SDF_CHANNEL_HDR(7).Description=[''];

SDF_CHANNEL_HDR(8).FieldIndex=8;
SDF_CHANNEL_HDR(8).BinaryIndex=[89 90];
SDF_CHANNEL_HDR(8).FieldName='weight';
SDF_CHANNEL_HDR(8).DataType='short';
SDF_CHANNEL_HDR(8).RangeUnits='0:3';
SDF_CHANNEL_HDR(8).Description=[''];
SDF_CHANNEL_HDR(8).Table=[...
{'0'} {'no weighting'};...
{'1'} {'A-weighting'};...
{'2'} {'B-weighting'};...
{'3'} {'C-weighting'};...
];

SDF_CHANNEL_HDR(30).FieldIndex=8;
SDF_CHANNEL_HDR(30).BinaryIndex=[91 94];
SDF_CHANNEL_HDR(30).FieldName='delayOld';
SDF_CHANNEL_HDR(30).DataType='float';
SDF_CHANNEL_HDR(30).RangeUnits='0:3';
SDF_CHANNEL_HDR(30).Description=['*** Prior to version 3.0'];

SDF_CHANNEL_HDR(9).FieldIndex=9;
SDF_CHANNEL_HDR(9).BinaryIndex=[93 96];
SDF_CHANNEL_HDR(9).FieldName='range';
SDF_CHANNEL_HDR(9).DataType='float';
SDF_CHANNEL_HDR(9).RangeUnits='unit is dBV, range is i.d. & includes overhead for scaling';
SDF_CHANNEL_HDR(9).Description=[''];

SDF_CHANNEL_HDR(10).FieldIndex=10;
SDF_CHANNEL_HDR(10).BinaryIndex=[99 100];
SDF_CHANNEL_HDR(10).FieldName='direction';
SDF_CHANNEL_HDR(10).DataType='short';
SDF_CHANNEL_HDR(10).RangeUnits='-9:9';
SDF_CHANNEL_HDR(10).Description=[''];
SDF_CHANNEL_HDR(10).Table=[...
 {'-9'} {'-TZ'};...
 {'-8'} {'-TY'};...
 {'-7'} {'-TX'};...
 {'-3'} {'-Z'};...
 {'-2'} {'-Y'};...
 {'-1'} {'-X'};...
 {'0'} {'no direction specified'};...
 {'1'} {'X'};...
 {'2'} {'Y'};...
 {'3'} {'Z'};...
 {'4'} {'R (radial)'};...
 {'5'} {'T (tangential � theta angle)'};...
 {'6'} {'P (tangential � phi angle)'};...
 {'7'} {'TX'};...
 {'8'} {'TY'};...
 {'9'} {'TZ'};...
];

SDF_CHANNEL_HDR(11).FieldIndex=11;
SDF_CHANNEL_HDR(11).BinaryIndex=[101 102];
SDF_CHANNEL_HDR(11).FieldName='pointNum';
SDF_CHANNEL_HDR(11).DataType='short';
SDF_CHANNEL_HDR(11).RangeUnits='-0:32676';
SDF_CHANNEL_HDR(11).Description=['test point on device under test'];

SDF_CHANNEL_HDR(12).FieldIndex=12;
SDF_CHANNEL_HDR(12).BinaryIndex=[103 104];
SDF_CHANNEL_HDR(12).FieldName='coupling';
SDF_CHANNEL_HDR(12).DataType='short';
SDF_CHANNEL_HDR(12).RangeUnits='0:1';
SDF_CHANNEL_HDR(12).Description=[''];
SDF_CHANNEL_HDR(12).Table=[...
 {'0'} {'DC'};...
 {'1'} {'AC'};...
];

SDF_CHANNEL_HDR(13).FieldIndex=13;
SDF_CHANNEL_HDR(13).BinaryIndex=[105 106];
SDF_CHANNEL_HDR(13).FieldName='overloaded';
SDF_CHANNEL_HDR(13).DataType='short';
SDF_CHANNEL_HDR(13).RangeUnits='0:1';
SDF_CHANNEL_HDR(13).Description=[''];
SDF_CHANNEL_HDR(13).Table=[...
 {'0'} {'no'};...
 {'1'} {'yes'};...
 {'2'} {'yes, error - undocumented value'};...
 {'3'} {'yes, error - undocumented value'};...
];

SDF_CHANNEL_HDR(14).FieldIndex=14;
SDF_CHANNEL_HDR(14).BinaryIndex=[107 116];
SDF_CHANNEL_HDR(14).FieldName='intLabel';
SDF_CHANNEL_HDR(14).DataType='char[10]';
SDF_CHANNEL_HDR(14).RangeUnits='i.d.';
SDF_CHANNEL_HDR(14).Description=...
['label for the instrument’s internal unit (such as V)'];

SDF_CHANNEL_HDR(15).FieldIndex=15;
SDF_CHANNEL_HDR(15).BinaryIndex=[117 138];
SDF_CHANNEL_HDR(15).FieldName='engUnit';
SDF_CHANNEL_HDR(15).DataType='struct';
SDF_CHANNEL_HDR(15).RangeUnits='see SDF_UNIT';
SDF_CHANNEL_HDR(15).Description=...
['engineering unit (EU) definition for this channel'];

SDF_CHANNEL_HDR(16).FieldIndex=16;
SDF_CHANNEL_HDR(16).BinaryIndex=[139 142];
SDF_CHANNEL_HDR(16).FieldName='int2engrUnit';
SDF_CHANNEL_HDR(16).DataType='float';
SDF_CHANNEL_HDR(16).RangeUnits='|10^34 (except 0)';
SDF_CHANNEL_HDR(16).Description=...
['EU correction factor. Divide internal-unit data ' ...
'by this value to get EU data'];

SDF_CHANNEL_HDR(17).FieldIndex=17;
SDF_CHANNEL_HDR(17).BinaryIndex=[143 146];
SDF_CHANNEL_HDR(17).FieldName='inputImpedance';
SDF_CHANNEL_HDR(17).DataType='float';
SDF_CHANNEL_HDR(17).RangeUnits='unit ohm, range i.d.';
SDF_CHANNEL_HDR(17).Description=['Input impedance'];

SDF_CHANNEL_HDR(18).FieldIndex=18;
SDF_CHANNEL_HDR(18).BinaryIndex=[147 148];
SDF_CHANNEL_HDR(18).FieldName='channelAttribute';
SDF_CHANNEL_HDR(18).DataType='short';
SDF_CHANNEL_HDR(18).RangeUnits='-99:3';
SDF_CHANNEL_HDR(18).Description=[''];
SDF_CHANNEL_HDR(18).Table=[...
 {'-99'} {'unknown attribute'};...
 {'0'} {'no attribute'};...
 {'1'} {'tach attribute'};...
 {'2'} {'reference attribute'};...
 {'3'} {'tach and reference attribute'};...
 {'4'} {'clockwise attribute'};...
];

SDF_CHANNEL_HDR(19).FieldIndex=19;
SDF_CHANNEL_HDR(19).BinaryIndex=[149 150];
SDF_CHANNEL_HDR(19).FieldName='aliasProtected';
SDF_CHANNEL_HDR(19).DataType='short';
SDF_CHANNEL_HDR(19).RangeUnits='0:1';
SDF_CHANNEL_HDR(19).Description=[''];
SDF_CHANNEL_HDR(19).Table=[...
 {'0'} {'data was not alias protected'};...
 {'1'} {'alias protected'};...
];

SDF_CHANNEL_HDR(20).FieldIndex=20;
SDF_CHANNEL_HDR(20).BinaryIndex=[151 152];
SDF_CHANNEL_HDR(20).FieldName='aliasProtected';
SDF_CHANNEL_HDR(20).DataType='short';
SDF_CHANNEL_HDR(20).RangeUnits='0:1';
SDF_CHANNEL_HDR(20).Description=[''];
SDF_CHANNEL_HDR(20).Table=[...
 {'0'} {'analog input channel'};...
 {'1'} {'digital input channel'};...
];

SDF_CHANNEL_HDR(21).FieldIndex=21;
SDF_CHANNEL_HDR(21).BinaryIndex=[153 160];
SDF_CHANNEL_HDR(21).FieldName='channelScale';
SDF_CHANNEL_HDR(21).DataType='double';
SDF_CHANNEL_HDR(21).RangeUnits='range is i.d.';
SDF_CHANNEL_HDR(21).Description=['see channelOffset below'];


SDF_CHANNEL_HDR(22).FieldIndex=22;
SDF_CHANNEL_HDR(22).BinaryIndex=[161 168];
SDF_CHANNEL_HDR(22).FieldName='channelOffset';
SDF_CHANNEL_HDR(22).DataType='double';
SDF_CHANNEL_HDR(22).RangeUnits='range is i.d.';
SDF_CHANNEL_HDR(22).Description=...
['when the data type is "short" or "long" the ' ...
'following formula will convert the data to volts:\n' ...
'Volts=channelOffset + ( channelScale ö Ydata)'];

SDF_CHANNEL_HDR(23).FieldIndex=23;
SDF_CHANNEL_HDR(23).BinaryIndex=[169 176];
SDF_CHANNEL_HDR(23).FieldName='gateBegin';
SDF_CHANNEL_HDR(23).DataType='double';
SDF_CHANNEL_HDR(23).RangeUnits='unit is sec, range is i.d.';
SDF_CHANNEL_HDR(23).Description=['Gated sweep start time'];

SDF_CHANNEL_HDR(24).FieldIndex=24;
SDF_CHANNEL_HDR(24).BinaryIndex=[177 184];
SDF_CHANNEL_HDR(24).FieldName='gateEnd';
SDF_CHANNEL_HDR(24).DataType='double';
SDF_CHANNEL_HDR(24).RangeUnits='unit is sec, range is i.d.';
SDF_CHANNEL_HDR(24).Description=['Gated sweep stop time'];

SDF_CHANNEL_HDR(25).FieldIndex=25;
SDF_CHANNEL_HDR(25).BinaryIndex=[185 192];
SDF_CHANNEL_HDR(25).FieldName='userDelay';
SDF_CHANNEL_HDR(25).DataType='double';
SDF_CHANNEL_HDR(25).RangeUnits='unit is sec, range is i.d.';
SDF_CHANNEL_HDR(25).Description=...
['User specified input channel time delay or line length (not trigger delay)'];

SDF_CHANNEL_HDR(26).FieldIndex=26;
SDF_CHANNEL_HDR(26).BinaryIndex=[193 200];
SDF_CHANNEL_HDR(26).FieldName='delay';
SDF_CHANNEL_HDR(26).DataType='double';
SDF_CHANNEL_HDR(26).RangeUnits='unit is sec, range is i.d.';
SDF_CHANNEL_HDR(26).Description=...
['amount of time between trigger event and start of data collection'];

SDF_CHANNEL_HDR(27).FieldIndex=27;
SDF_CHANNEL_HDR(27).BinaryIndex=[201 208];
SDF_CHANNEL_HDR(27).FieldName='carrierFreq';
SDF_CHANNEL_HDR(27).DataType='double';
SDF_CHANNEL_HDR(27).RangeUnits='unit is Hz, range is i.d.';
SDF_CHANNEL_HDR(27).Description=['carrier frequency for demodulated data'];

SDF_CHANNEL_HDR(28).FieldIndex=28;
SDF_CHANNEL_HDR(28).BinaryIndex=[209 210];
SDF_CHANNEL_HDR(28).FieldName='channelNumber';
SDF_CHANNEL_HDR(28).DataType='short';
SDF_CHANNEL_HDR(28).RangeUnits='0:32767';
SDF_CHANNEL_HDR(28).Description=['zero-based channel number'];

SDF_CHANNEL_HDR(29).FieldIndex=29;
SDF_CHANNEL_HDR(29).BinaryIndex=[211 212];
SDF_CHANNEL_HDR(29).FieldName='channelModule';
SDF_CHANNEL_HDR(29).DataType='short';
SDF_CHANNEL_HDR(29).RangeUnits='0:32767';
SDF_CHANNEL_HDR(29).Description=['zero-based channel module'];
end

function  SDF_COMMENT_HDR =  SDF_COMMENT_HDR_Template()
%SDF_COMMENT_HDR_TEMPLATE - Provides a empty template containing the requested header information
%
%                 - SDF_COMMENT_HDR: Contains text information associated with the file.
%
% Syntax: [retval] = SDF_COMMENT_HDR_TEMPLATE()
%
% Inputs: none
%
% Outputs:
%   retval - Extracted SDF header information.
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: N/A
% © Copyright 2016 Commonwealth of Australia, represented by the Department of Defence

SDF_COMMENT_HDR(1).FieldIndex=1;
SDF_COMMENT_HDR(1).BinaryIndex=[1 2];
SDF_COMMENT_HDR(1).FieldName='recordType';
SDF_COMMENT_HDR(1).DataType = 'short';
SDF_COMMENT_HDR(1).RangeUnits='20';

SDF_COMMENT_HDR(2).FieldIndex=2;
SDF_COMMENT_HDR(2).BinaryIndex=[3 6];
SDF_COMMENT_HDR(2).FieldName='recordSize';
SDF_COMMENT_HDR(2).DataType = 'long';
SDF_COMMENT_HDR(2).RangeUnits='variable';

SDF_COMMENT_HDR(3).FieldIndex=3;
SDF_COMMENT_HDR(3).BinaryIndex=[7 10];
SDF_COMMENT_HDR(3).FieldName='unique_record';
SDF_COMMENT_HDR(3).DataType = 'long';
SDF_COMMENT_HDR(3).RangeUnits='-1:(231)-1';
SDF_COMMENT_HDR(3).Description=...
['byte offset from the beginning of the file to a ' ...
'record containing an instrument-specific comment header. May ' ...
'be ignored if recalled into a different type instrument.'];

SDF_COMMENT_HDR(4).FieldIndex=4;
SDF_COMMENT_HDR(4).BinaryIndex=[11 14];
SDF_COMMENT_HDR(4).FieldName='headersize';
SDF_COMMENT_HDR(4).DataType = 'long';
SDF_COMMENT_HDR(4).RangeUnits='24';
SDF_COMMENT_HDR(4).Description=...
['size of the header portion of this record (excluding the comment text).'];

SDF_COMMENT_HDR(5).FieldIndex=5;
SDF_COMMENT_HDR(5).BinaryIndex=[15 18];
SDF_COMMENT_HDR(5).FieldName='comment_bytes';
SDF_COMMENT_HDR(5).DataType = 'long';
SDF_COMMENT_HDR(5).RangeUnits='-1:recordSize-headerSize';
SDF_COMMENT_HDR(5).Description=...
['size of comment (in bytes). This size may be ' ...
'smaller than the comment text area. If the size of the text is -1, ' ...
'then the entire comment text area is valid (or until an end-of-text ' ...
'marker is found).'];

SDF_COMMENT_HDR(6).FieldIndex=6;
SDF_COMMENT_HDR(6).BinaryIndex=[19 20];
SDF_COMMENT_HDR(6).FieldName='comment_type';
SDF_COMMENT_HDR(6).DataType = 'short';
SDF_COMMENT_HDR(6).RangeUnits='0:0';
SDF_COMMENT_HDR(6).Description=...
['type of comment data'];

SDF_COMMENT_HDR(7).FieldIndex=7;
SDF_COMMENT_HDR(7).BinaryIndex=[21 22];
SDF_COMMENT_HDR(7).FieldName='scope_type';
SDF_COMMENT_HDR(7).DataType = 'short';
SDF_COMMENT_HDR(7).RangeUnits='0:4';
SDF_COMMENT_HDR(7).Description=...
['tells which type of header the comment applies to'];
SDF_COMMENT_HDR(7).Table=[...
{'0'} {'entire file'};...
{'1'} {'SDF_DATA_HDR'};...
{'2'} {'SDF_VECTOR_HDR'};...
{'3'} {'SDF_CHANNEL_HDR'};...
{'4'} {'SDF_SCAN_STRUCT'};...
];

SDF_COMMENT_HDR(8).FieldIndex=8;
SDF_COMMENT_HDR(8).BinaryIndex=[23 24];
SDF_COMMENT_HDR(8).FieldName='scope_info';
SDF_COMMENT_HDR(8).DataType = 'short';
SDF_COMMENT_HDR(8).RangeUnits='-1:32767';
SDF_COMMENT_HDR(8).Description=...
['the index of the header associated with the ' ...
'scope_type (-1 = no specific header)'];
end

function  SDF_DATA_HDR=SDF_DATA_HDR_Template()
%SDF_DATA_HDR_TEMPLATE - Provides a empty template containing the requested header information
%
%                 - SDF_DATA_HDR:    Tells you how reconstruct one block of measurement results 
%                                    (x- and y-axis values for every point of every trace).
%
% Syntax: [retval] = SDF_DATA_HDR_TEMPLATE()
%
% Inputs: none
%
% Outputs:
%   retval - Extracted SDF header information.
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: N/A
% © Copyright 2016 Commonwealth of Australia, represented by the Department of Defence

SDF_DATA_HDR(1).FieldIndex=1;
SDF_DATA_HDR(1).BinaryIndex=[1 2];
SDF_DATA_HDR(1).FieldName='recordType';
SDF_DATA_HDR(1).DataType='short';
SDF_DATA_HDR(1).RangeUnits='12';

SDF_DATA_HDR(2).FieldIndex=2;
SDF_DATA_HDR(2).BinaryIndex=[3 6];
SDF_DATA_HDR(2).FieldName='recordSize';
SDF_DATA_HDR(2).DataType='long';
SDF_DATA_HDR(2).RangeUnits='148 bytes';

SDF_DATA_HDR(3).FieldIndex=3;
SDF_DATA_HDR(3).BinaryIndex=[7 10];
SDF_DATA_HDR(3).FieldName='unique_record';
SDF_DATA_HDR(3).DataType='long';
SDF_DATA_HDR(3).RangeUnits='-1:2^31-1';
SDF_MEAS_HDR(3).Description=['byte offset from the beginning of the file ' ...
'to a record containing an instrument-specific data header. ' ...
'May be ignored if recalled into a different type instrument.'];

SDF_DATA_HDR(4).FieldIndex=4;
SDF_DATA_HDR(4).BinaryIndex=[11 26];
SDF_DATA_HDR(4).FieldName='dataTitle';
SDF_DATA_HDR(4).DataType='char[16]';
SDF_DATA_HDR(4).RangeUnits='i.d.';
SDF_MEAS_HDR(4).Description=['instrument- or user-supplied name ' ...
'for data type.'];

SDF_DATA_HDR(5).FieldIndex=5;
SDF_DATA_HDR(5).BinaryIndex=[27 28];
SDF_DATA_HDR(5).FieldName='domain';
SDF_DATA_HDR(5).DataType='short';
SDF_DATA_HDR(5).RangeUnits='-99, 0:6';
SDF_MEAS_HDR(5).Description=[''];
SDF_DATA_HDR(5).Table=[...
 {'-99'} {'unknown'};...
 {'0'} {'frequency'};...
 {'1'} {'time'};...
 {'2'} {'amplitude'};...
 {'3'} {'RPM'};...
 {'4'} {'order'};...
 {'5'} {'channel'};...
 {'6'} {'octave'};...
];

SDF_DATA_HDR(6).FieldIndex=6;
SDF_DATA_HDR(6).BinaryIndex=[29 30];
SDF_DATA_HDR(6).FieldName='dataType';
SDF_DATA_HDR(6).DataType='short';
SDF_DATA_HDR(6).RangeUnits='-99, 0:76';
SDF_MEAS_HDR(6).Description=[''];
SDF_DATA_HDR(6).Table=[...
 {'-99'} {'unknown'};...
 {'0'} {'time'};...
 {'1'} {'linear spectrum'};...
 {'2'} {'auto-power spectrum'};...
 {'3'} {'cross-power spectrum'};...
 {'4'} {'frequency response'};...
 {'5'} {'auto-correlation'};...
 {'6'} {'cross-correlation'};...
 {'7'} {'impulse response'};...
 {'8'} {'ordinary coherence'};...
 {'9'} {'partial coherence'};...
 {'10'} {'multiple coherence'};...
 {'11'} {'full octave'};...
 {'12'} {'third octave'};...
 {'13'} {'convolution'};...
 {'14'} {'histogram'};...
 {'15'} {'probability density function'};...
 {'16'} {'cumulative density function'};...
 {'17'} {'power spectrum order tracking'};...
 {'18'} {'composite power tracking'};...
 {'19'} {'phase order tracking'};...
 {'20'} {'rpm spectral'};...
 {'21'} {'order ratio'};...
 {'22'} {'orbit'};...
 {'23'} {'HP 35650 series calibration'};...
 {'24'} {'sine rms pwr data'};...
 {'25'} {'sine variance data'};...
 {'26'} {'sine range data'};...
 {'27'} {'sine settle time data'};...
 {'28'} {'sine integ time data'};...
 {'29'} {'sine source data'};...
 {'30'} {'sine overload data'};...
 {'31'} {'sine linear data'};...
 {'32'} {'synthesis'};...
 {'33'} {'curve fit weighting function'};...
 {'34'} {'frequency corrections (for capture)'};...
 {'35'} {'all pass time data'};...
 {'36'} {'norm reference data'};...
 {'37'} {'tachometer data'};...
 {'38'} {'limit line data'};...
 {'39'} {'twelfth octave data'};...
 {'40'} {'S11 data'};...
 {'41'} {'S21 data'};...
 {'42'} {'S12 data'};...
 {'43'} {'S22 data'};...
 {'44'} {'PSD data'};...
 {'45'} {'decimated time data'};...
 {'46'} {'overload data'};...
 {'47'} {'compressed time data'};...
 {'48'} {'external trigger data'};...
 {'49'} {'pressure data'};...
 {'50'} {'intensity data'};...
 {'51'} {'PI index data'};...
 {'52'} {'velocity data'};...
 {'53'} {'PV index data'};...
 {'54'} {'sound power data'};...
 {'55'} {'field indicator data'};...
 {'56'} {'partial power data'};...
 {'57'} {'Ln 1 data'};...
 {'58'} {'Ln 10 data'};...
 {'59'} {'Ln 50 data'};...
 {'60'} {'Ln 90 data'};...
 {'61'} {'Ln 99 data'};...
 {'62'} {'Ln user data'};...
 {'63'} {'T20 data'};...
 {'64'} {'T30 data'};...
 {'65'} {'RT60 data'};...
 {'66'} {'average count data'};...
 {'68'} {'IQ measured time'};...
 {'69'} {'IQ measured spectrum'};...
 {'70'} {'IQ reference time'};...
 {'71'} {'IQ reference spectrum'};...
 {'72'} {'IQ error magnitude'};...
 {'73'} {'IQ error phase'};...
 {'74'} {'IQ error vector time'};...
 {'75'} {'IQ error vector spectrum'};...
 {'76'} {'symbol table data'};...
];

SDF_DATA_HDR(7).FieldIndex=7;
SDF_DATA_HDR(7).BinaryIndex=[31 32];
SDF_DATA_HDR(7).FieldName='num_of_pointsOld';
SDF_DATA_HDR(7).DataType='short';
SDF_DATA_HDR(7).RangeUnits='';
SDF_MEAS_HDR(7).Description=['*** Prior to version 3.0'];

SDF_DATA_HDR(8).FieldIndex=8;
SDF_DATA_HDR(8).BinaryIndex=[33 34];
SDF_DATA_HDR(8).FieldName='last_valid_indexOld';
SDF_DATA_HDR(8).DataType='short';
SDF_DATA_HDR(8).RangeUnits='';
SDF_MEAS_HDR(8).Description=['*** Prior to version 3.0'];

SDF_DATA_HDR(9).FieldIndex=9;
SDF_DATA_HDR(9).BinaryIndex=[35 38];
SDF_DATA_HDR(9).FieldName='abscissa_firstXOld';
SDF_DATA_HDR(9).DataType='float';
SDF_DATA_HDR(9).RangeUnits='';
SDF_MEAS_HDR(9).Description=['** Prior to version 2.0'];

SDF_DATA_HDR(10).FieldIndex=10;
SDF_DATA_HDR(10).BinaryIndex=[39 42];
SDF_DATA_HDR(10).FieldName='abscissa_deltaXOld';
SDF_DATA_HDR(10).DataType='float';
SDF_DATA_HDR(10).RangeUnits='';
SDF_MEAS_HDR(10).Description=['** Prior to version 2.0'];

SDF_DATA_HDR(11).FieldIndex=11;
SDF_DATA_HDR(11).BinaryIndex=[43 44];
SDF_DATA_HDR(11).FieldName='xResolution_type';
SDF_DATA_HDR(11).DataType='short';
SDF_DATA_HDR(11).RangeUnits='0:4';
SDF_MEAS_HDR(11).Description=['tells you how to find x-axis values for this ' ...
'Data Header record’s traces.\n' ...
'0=linear —calculate values from abscissa_firstX and abscissa_deltaX\n' ...
'1=logarithmic—calculate values from abscissa_firstX and abscissa_deltaX\n' ...
'2=arbitrary, one per file—find values in the X-axis Data record; same ' ...
'vector used for every trace in the measurement file\n' ...
'3=arbitrary, one per data type—find values in the X-axis Data record; ' ...
'same x-axis vector used for each trace associated with this record\n' ...
'4=arbitrary, one per trace—find values in the X-axis Data record; ' ...
'unique x-axis vector for each trace associated with this record'];
SDF_DATA_HDR(11).Table=[...
 {'0'} {'linear'};...
 {'1'} {'logarithmic'};...
 {'2'} {'arbitrary, one per file'};...
 {'3'} {'arbitrary, one per data type'};...
 {'4'} {'arbitrary, one per trace'};...
];

SDF_DATA_HDR(12).FieldIndex=12;
SDF_DATA_HDR(12).BinaryIndex=[45 46];
SDF_DATA_HDR(12).FieldName='xdata_type';
SDF_DATA_HDR(12).DataType='short';
SDF_DATA_HDR(12).RangeUnits='1:4';
SDF_MEAS_HDR(12).Description=...
['tells you the size and format of each x-axis value.' ...
'1=short (two-byte, binary-encoded integer)\n' ...
'2=long (four-byte, binary-encoded integer)\n' ...
'3=float (four-byte, binary floating-point number)\n' ...
'4=double (eight-byte, binary floating-point number)\n' ...
'This field is only valid if xResolution_type is 2, 3, or 4.'];
SDF_DATA_HDR(12).Table=[...
 {'1'} {'short'};...
 {'2'} {'long'};...
 {'3'} {'float'};...
 {'4'} {'double'};...
];

SDF_DATA_HDR(13).FieldIndex=13;
SDF_DATA_HDR(13).BinaryIndex=[47 48];
SDF_DATA_HDR(13).FieldName='xPerPoint';
SDF_DATA_HDR(13).DataType='short';
SDF_DATA_HDR(13).RangeUnits='0:32767';
SDF_MEAS_HDR(13).Description=...
['number of x-axis values per each trace point. ' ...
'This field is only valid if xResolution_type is 2, 3, or 4.'];

SDF_DATA_HDR(14).FieldIndex=14;
SDF_DATA_HDR(14).BinaryIndex=[49 50];
SDF_DATA_HDR(14).FieldName='ydata_type';
SDF_DATA_HDR(14).DataType='short';
SDF_DATA_HDR(14).RangeUnits='1:4';
SDF_MEAS_HDR(14).Description=...
['tells you the size and format of each y-axis value.\n' ...
'NOTE: If yIsComplex=1, both the real and imaginary ' ...
'components of each y value require the number of bytes specified here.\n' ...
'1=short (two-byte, binary-encoded integer)\n' ...
'2=long (four-byte, binary-encoded integer)\n' ...
'3=float (four-byte, binary floating-point number)\n' ...
'4=double (eight-byte, binary floating-point number)\n'];
SDF_DATA_HDR(14).Table=[...
 {'1'} {'short'};...
 {'2'} {'long'};...
 {'3'} {'float'};...
 {'4'} {'double'};...
];

SDF_DATA_HDR(15).FieldIndex=15;
SDF_DATA_HDR(15).BinaryIndex=[51 52];
SDF_DATA_HDR(15).FieldName='yPerPoint';
SDF_DATA_HDR(15).DataType='short';
SDF_DATA_HDR(15).RangeUnits='0:32767';
SDF_MEAS_HDR(15).Description=...
['number of y-axis values per each trace point.\n' ...
'NOTE: A value containing both real and imaginary ' ...
'components is still considered a single value.'];

SDF_DATA_HDR(16).FieldIndex=16;
SDF_DATA_HDR(16).BinaryIndex=[53 54];
SDF_DATA_HDR(16).FieldName='yIsComplex';
SDF_DATA_HDR(16).DataType='short';
SDF_DATA_HDR(16).RangeUnits='0:1';
SDF_MEAS_HDR(16).Description=[''];
SDF_DATA_HDR(16).Table=[...
 {'0'} {'each y value has only a real component'};...
 {'1'} {'each y value has both a real and an imaginary component'};...
];

SDF_DATA_HDR(17).FieldIndex=17;
SDF_DATA_HDR(17).BinaryIndex=[55 56];
SDF_DATA_HDR(17).FieldName='yIsNormalised';
SDF_DATA_HDR(17).DataType='short';
SDF_DATA_HDR(17).RangeUnits='0:1';
SDF_MEAS_HDR(17).Description=['0=not normalized\n' ...
'1=normalized (all y values fall between 0.0 and 1.0 and ' ...
'are unitless, for example coherence or power spectrum).'];
SDF_DATA_HDR(17).Table=[...
 {'0'} {'not normalised'};...
 {'1'} {'normalized'};...
];

SDF_DATA_HDR(18).FieldIndex=18;
SDF_DATA_HDR(18).BinaryIndex=[57 58];
SDF_DATA_HDR(18).FieldName='yIsPowerData';
SDF_DATA_HDR(18).DataType='short';
SDF_DATA_HDR(18).RangeUnits='0:1';
SDF_MEAS_HDR(18).Description=...
['0=not power data (for example, linear spectrum)\n' ...
'1=power data (for example, auto-power spectrum)'];
SDF_DATA_HDR(18).Table=[...
 {'0'} {'not power data'};...
 {'1'} {'power data'};...
];

SDF_DATA_HDR(19).FieldIndex=19;
SDF_DATA_HDR(19).BinaryIndex=[59 60];
SDF_DATA_HDR(19).FieldName='yIsValid';
SDF_DATA_HDR(19).DataType='short';
SDF_DATA_HDR(19).RangeUnits='0:1';
SDF_MEAS_HDR(19).Description=[''];
SDF_DATA_HDR(19).Table=[...
{'0'} {'non valid'};...
{'1'} {'valid'};...
];

SDF_DATA_HDR(20).FieldIndex=20;
SDF_DATA_HDR(20).BinaryIndex=[61 64];
SDF_DATA_HDR(20).FieldName='first_VECTOR_recordNum';
SDF_DATA_HDR(20).DataType='long';
SDF_DATA_HDR(20).RangeUnits='0:(num_of_VECTOR_record) − 1(SDF_FILE_HDR)';
SDF_MEAS_HDR(20).Description=['first Vector Header record ' ...
'belonging to this Data Header record'];

SDF_DATA_HDR(21).FieldIndex=21;
SDF_DATA_HDR(21).BinaryIndex=[65 66];
SDF_DATA_HDR(21).FieldName='total_rows';
SDF_DATA_HDR(21).DataType='short';
SDF_DATA_HDR(21).RangeUnits='1:32767';
SDF_MEAS_HDR(21).Description=['used to determine the number of traces ' ...
'associated with this Data header record; just multiply ' ...
'total_rows by total_columns'];

SDF_DATA_HDR(22).FieldIndex=22;
SDF_DATA_HDR(22).BinaryIndex=[67 68];
SDF_DATA_HDR(22).FieldName='total_cols';
SDF_DATA_HDR(22).DataType='short';
SDF_DATA_HDR(22).RangeUnits='1:32767';
SDF_MEAS_HDR(22).Description=['used to determine the number of traces ' ...
'associated with this Data header record; just multiply ' ...
'total_rows by total_columns'];

SDF_DATA_HDR(23).FieldIndex=23;
SDF_DATA_HDR(23).BinaryIndex=[69 90];
SDF_DATA_HDR(23).FieldName='xUnit';
SDF_DATA_HDR(23).DataType='struct';
SDF_DATA_HDR(23).RangeUnits='see SDF_UNIT';
SDF_MEAS_HDR(23).Description=['engineering unit used for x-axis values.'];

SDF_DATA_HDR(24).FieldIndex=24;
SDF_DATA_HDR(24).BinaryIndex=[91 92];
SDF_DATA_HDR(24).FieldName='yUnitValid';
SDF_DATA_HDR(24).DataType='short';
SDF_DATA_HDR(24).RangeUnits='0:1';
SDF_MEAS_HDR(24).Description=['yUnit field in this record is valid'];
SDF_DATA_HDR(24).Table=[...
{'0'} {'use Channel Header record’s engUnit field for y-axis units.'};...
{'1'} {'use this record’s yUnit field for y-axis unit.'};...
];

SDF_DATA_HDR(25).FieldIndex=25;
SDF_DATA_HDR(25).BinaryIndex=[93 114];
SDF_DATA_HDR(25).FieldName='yUnit';
SDF_DATA_HDR(25).DataType='struct';
SDF_DATA_HDR(25).RangeUnits='see SDF_UNIT';
SDF_MEAS_HDR(25).Description=['engineering unit used for y-axis values.'];

SDF_DATA_HDR(26).FieldIndex=26;
SDF_DATA_HDR(26).BinaryIndex=[115 122];
SDF_DATA_HDR(26).FieldName='abscissa_firstX';
SDF_DATA_HDR(26).DataType='double';
SDF_DATA_HDR(26).RangeUnits='|10^34';
SDF_MEAS_HDR(26).Description=...
['firstX—x-axis value of first point. This field is only ' ...
'valid if xResolution_type is 0 or 1.'];

SDF_DATA_HDR(27).FieldIndex=27;
SDF_DATA_HDR(27).BinaryIndex=[123 130];
SDF_DATA_HDR(27).FieldName='abscissa_deltaX';
SDF_DATA_HDR(27).DataType='double';
SDF_DATA_HDR(27).RangeUnits='|10^34';
SDF_MEAS_HDR(27).Description=['spacing between x-axis points.\n' ...
'xn=x(n-1) + abscissa_deltaX (xResolution_type 0)\n' ...
'xn=x(n-1)öabscissa_deltaX (xResolution_type 1)\n' ...
'This field is only valid if xResolution_type is 0 or 1.'];

SDF_DATA_HDR(28).FieldIndex=28;
SDF_DATA_HDR(28).BinaryIndex=[131 132];
SDF_DATA_HDR(28).FieldName='scanData';
SDF_DATA_HDR(28).DataType='short';
SDF_DATA_HDR(28).RangeUnits='0:1';
SDF_MEAS_HDR(28).Description=['indicates whether the data is scanned'];
SDF_DATA_HDR(28).Table=[...
{'0'} {'This data header is associated with non-scanned data'};...
{'1'} {'This data header is associated with the scan structure'};...
];

SDF_DATA_HDR(29).FieldIndex=29;
SDF_DATA_HDR(29).BinaryIndex=[133 134];
SDF_DATA_HDR(29).FieldName='windowApplied';
SDF_DATA_HDR(29).DataType='short';
SDF_DATA_HDR(29).RangeUnits='0:1';
SDF_MEAS_HDR(29).Description=['indicates whether the windows indicated ' ...
'have already been applied to the data'];
SDF_DATA_HDR(29).Table=[...
{'0'} {'windows have not been applied'};...
{'1'} {'windows have been applied'};...
];

SDF_DATA_HDR(30).FieldIndex=30;
SDF_DATA_HDR(30).BinaryIndex=[135 138];
SDF_DATA_HDR(30).FieldName='num_of_points';
SDF_DATA_HDR(30).DataType='long';
SDF_DATA_HDR(30).RangeUnits='0:(2^31)-1';
SDF_MEAS_HDR(30).Description=['number of discrete points in each trace ' ...
'associated with this record.'];

SDF_DATA_HDR(31).FieldIndex=31;
SDF_DATA_HDR(31).BinaryIndex=[139 142];
SDF_DATA_HDR(31).FieldName='last_valid_index';
SDF_DATA_HDR(31).DataType='long';
SDF_DATA_HDR(31).RangeUnits='0:(num_of_points)-1';
SDF_MEAS_HDR(31).Description=['last point containing valid data.'];

SDF_DATA_HDR(32).FieldIndex=32;
SDF_DATA_HDR(32).BinaryIndex=[143 144];
SDF_DATA_HDR(32).FieldName='overSampleFactor';
SDF_DATA_HDR(32).DataType='short';
SDF_DATA_HDR(32).RangeUnits='0:32767';
SDF_MEAS_HDR(32).Description=['Usually 1\n' ...
'> 1=the data has been low-pass filtered but not decimated.'];

SDF_DATA_HDR(33).FieldIndex=33;
SDF_DATA_HDR(33).BinaryIndex=[145 146];
SDF_DATA_HDR(33).FieldName='multiPassMode';
SDF_DATA_HDR(33).DataType='short';
SDF_DATA_HDR(33).RangeUnits='0:5';
SDF_MEAS_HDR(33).Description=...
['"Multi-pass" refers to a mode where data for ' ...
'multiple frequency spans is interleaved.'];
SDF_DATA_HDR(33).Table=[...
{'0'} {'non multipass data'};...
{'1'} {'multi-pass, corresponding the HP 3565 gate array modes'};...
{'2'} {'multi-pass, corresponding the HP 3565 gate array modes'};...
{'3'} {'multi-pass, corresponding the HP 3565 gate array modes'};...
{'4'} {'multi-pass, corresponding the HP 3565 gate array modes'};...
{'5'} {'future multi-pass mode'};...
];

SDF_DATA_HDR(34).FieldIndex=34;
SDF_DATA_HDR(34).BinaryIndex=[147 148];
SDF_DATA_HDR(34).FieldName='multiPassDecimations';
SDF_DATA_HDR(34).DataType='short';
SDF_DATA_HDR(34).RangeUnits='0:32767';
SDF_MEAS_HDR(34).Description=...
['> 0=the number of decimations included in the multi-pass data.  '];
end

function SDF_FILE_HDR=SDF_FILE_HDR_Template()
%SDF_FILE_HDR_TEMPLATE - Provides a empty template containing the requested header information
% 
%                 - SDF_FILE_HDR:    Provides an index to the file.
%
% Syntax: [retval] = SDF_FILE_HDR_TEMPLATE()
%
% Inputs: none
%
% Outputs:
%   retval - Extracted SDF header information.
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: N/A
% © Copyright 2016 Commonwealth of Australia, represented by the Department of Defence

SDF_FILE_HDR(1).FieldIndex=1;
SDF_FILE_HDR(1).BinaryIndex=[1 2];
SDF_FILE_HDR(1).FieldName=['recordType'];
SDF_FILE_HDR(1).DataType='short';
SDF_FILE_HDR(1).RangeUnits='10';

SDF_FILE_HDR(2).FieldIndex=2;
SDF_FILE_HDR(2).BinaryIndex=[3 6];
SDF_FILE_HDR(2).FieldName=['recordSize'];
SDF_FILE_HDR(2).DataType='long';
SDF_FILE_HDR(2).RangeUnits='80 bytes';

SDF_FILE_HDR(3).FieldIndex=3;
SDF_FILE_HDR(3).BinaryIndex=[7 8];
SDF_FILE_HDR(3).FieldName=['revisionNum'];
SDF_FILE_HDR(3).DataType='short';
SDF_FILE_HDR(3).RangeUnits='0:32767';
SDF_FILE_HDR(3).Description='measurement file version number.';

SDF_FILE_HDR(4).FieldIndex=4;
SDF_FILE_HDR(4).BinaryIndex=[9 10];
SDF_FILE_HDR(4).FieldName=['applic'];
SDF_FILE_HDR(4).DataType='short';
SDF_FILE_HDR(4).RangeUnits='-99:32767';
SDF_FILE_HDR(4).Description='file saved from this instrument or application.';
SDF_FILE_HDR(4).Table=...
    [{'-1'} {'HP VISTA'};...
    {'-2'} {'HP SINE'};...
    {'-3'} {'HP 35660A'};...
    {'-4'} {'HP 3562A, HP 3563A'};...
    {'-5'} {'HP 3588A'};...
    {'-6'} {'HP 3589A'};...
    {'-99'} {'unknown'};...
    {'1'} {'HP 3566A, HP 3567A'};...
    {'2'} {'HP 35665A'};...
    {'3'} {'HP 3560A'};...
    {'4'} {'HP 89410A, HP 89440A'};...
    {'7'} {'HP 35635R'};...
    {'8'} {'HP 35654A-S1A'};...
    {'9'} {'HP 3569A'};...
    {'10'} {'HP 35670A'};...
    {'11'} {'HP 3587S'};...
    ];

SDF_FILE_HDR(5).FieldIndex=5;
SDF_FILE_HDR(5).BinaryIndex=[11 12];
SDF_FILE_HDR(5).FieldName=['yearStamp'];
SDF_FILE_HDR(5).DataType='short';
SDF_FILE_HDR(5).RangeUnits='0:9999';
SDF_FILE_HDR(5).Description = 'year at measurement start.';


SDF_FILE_HDR(6).FieldIndex=6;
SDF_FILE_HDR(6).BinaryIndex=[13 14];
SDF_FILE_HDR(6).FieldName=['monthDayStamp'];
SDF_FILE_HDR(6).DataType='short';
SDF_FILE_HDR(6).RangeUnits='0:1231';
SDF_FILE_HDR(6).Description = ['month at measurement start.\n' ...
'Encoded as (month ö 100) + day\n' ...
'For example, November 9 = 1109'];

SDF_FILE_HDR(7).FieldIndex=7;
SDF_FILE_HDR(7).BinaryIndex=[15 16];
SDF_FILE_HDR(7).FieldName=['hourMinStamp'];
SDF_FILE_HDR(7).DataType='short';
SDF_FILE_HDR(7).RangeUnits='0:2359';
SDF_FILE_HDR(7).Description = ['time at measurement start.' ...
'Encoded as (hour ö 100) + minute' ...
'For example, 14:45 = 1445'];

SDF_FILE_HDR(8).FieldIndex=8;
SDF_FILE_HDR(8).BinaryIndex=[17 24];
SDF_FILE_HDR(8).FieldName=['applicVer'];
SDF_FILE_HDR(8).DataType='char[8]';
SDF_FILE_HDR(8).RangeUnits='i.d.';
SDF_FILE_HDR(8).Description = 'software or firmware version number.';

SDF_FILE_HDR(9).FieldIndex=9;
SDF_FILE_HDR(9).BinaryIndex=[25 26];
SDF_FILE_HDR(9).FieldName=['num_of_DATA_HDR_record'];
SDF_FILE_HDR(9).DataType='short';
SDF_FILE_HDR(9).RangeUnits='1:32767';
SDF_FILE_HDR(9).Description = 'total Data Header records.';

SDF_FILE_HDR(10).FieldIndex=10;
SDF_FILE_HDR(10).BinaryIndex=[27 28];
SDF_FILE_HDR(10).FieldName=['num_of_VECTOR_record'];
SDF_FILE_HDR(10).DataType='short';
SDF_FILE_HDR(10).RangeUnits='1:32767';
SDF_FILE_HDR(10).Description = 'total Vector Header records.';

SDF_FILE_HDR(11).FieldIndex=11;
SDF_FILE_HDR(11).BinaryIndex=[29 30];
SDF_FILE_HDR(11).FieldName=['num_of_CHANNEL_record'];
SDF_FILE_HDR(11).DataType='short';
SDF_FILE_HDR(11).RangeUnits='1:32767';
SDF_FILE_HDR(11).Description = 'total Channel Header records.';

SDF_FILE_HDR(12).FieldIndex=12;
SDF_FILE_HDR(12).BinaryIndex=[31 32];
SDF_FILE_HDR(12).FieldName=['num_of_UNIQUE_record'];
SDF_FILE_HDR(12).DataType='short';
SDF_FILE_HDR(12).RangeUnits='0:32767';
SDF_FILE_HDR(12).Description = 'total Unique records.';

SDF_FILE_HDR(13).FieldIndex=13;
SDF_FILE_HDR(13).BinaryIndex=[33 34];
SDF_FILE_HDR(13).FieldName=['num_of_SCAN_STRUCT_record'];
SDF_FILE_HDR(13).DataType='short';
SDF_FILE_HDR(13).RangeUnits='0:1';
SDF_FILE_HDR(13).Description = 'total Scan Structure records.';

SDF_FILE_HDR(14).FieldIndex=14;
SDF_FILE_HDR(14).BinaryIndex=[35 36];
SDF_FILE_HDR(14).FieldName=['num_of_XDATA_record'];
SDF_FILE_HDR(14).DataType='short';
SDF_FILE_HDR(14).RangeUnits='0:1';
SDF_FILE_HDR(14).Description = 'total X-axis Data records.';

SDF_FILE_HDR(15).FieldIndex=15;
SDF_FILE_HDR(15).BinaryIndex=[37 40];
SDF_FILE_HDR(15).FieldName=['offset_of_DATA_HDR_record'];
SDF_FILE_HDR(15).DataType='long';
SDF_FILE_HDR(15).RangeUnits='-1:(2^31)-1';
SDF_FILE_HDR(15).Description = ['first Data Header record’s byte offset ' ...
'from beginning of file.'];
 
SDF_FILE_HDR(16).FieldIndex=16;
SDF_FILE_HDR(16).BinaryIndex=[41 44];
SDF_FILE_HDR(16).FieldName=['offset_of_VECTOR_record'];
SDF_FILE_HDR(16).DataType='long';
SDF_FILE_HDR(16).RangeUnits='-1:(2^31)-1';
SDF_FILE_HDR(16).Description = ['first Vector Header record’s byte offset ' ...
'from beginning of file.'];

SDF_FILE_HDR(17).FieldIndex=17;
SDF_FILE_HDR(17).BinaryIndex=[45 48];
SDF_FILE_HDR(17).FieldName=['offset_of_CHANNEL_record'];
SDF_FILE_HDR(17).DataType='long';
SDF_FILE_HDR(17).RangeUnits='-1:(2^31)-1';
SDF_FILE_HDR(17).Description = ['first Channel Header record’s byte offset ' ...
'from beginning of file.'];

SDF_FILE_HDR(18).FieldIndex=18;
SDF_FILE_HDR(18).BinaryIndex=[49 52];
SDF_FILE_HDR(18).FieldName=['offset_of_UNIQUE_record'];
SDF_FILE_HDR(18).DataType='long';
SDF_FILE_HDR(18).RangeUnits='-1:(2^31)-1';
SDF_FILE_HDR(18).Description = ['first Unique record’s byte offset ' ...
'from beginning of file.'];

SDF_FILE_HDR(19).FieldIndex=19;
SDF_FILE_HDR(19).BinaryIndex=[53 56];
SDF_FILE_HDR(19).FieldName=['offset_of_SCAN_STRUCT_record'];
SDF_FILE_HDR(19).DataType='long';
SDF_FILE_HDR(19).RangeUnits='-1:(2^31)-1';
SDF_FILE_HDR(19).Description = ['Scan Structure record’s byte offset ' ...
'from beginning of file.'];

SDF_FILE_HDR(20).FieldIndex=20;
SDF_FILE_HDR(20).BinaryIndex=[57 60];
SDF_FILE_HDR(20).FieldName=['offset_of_XDATA_record'];
SDF_FILE_HDR(20).DataType='long';
SDF_FILE_HDR(20).RangeUnits='-1:(2^31)-1';
SDF_FILE_HDR(20).Description = ['X-axis Data record’s byte offset ' ...
'from beginning of file.'];

SDF_FILE_HDR(21).FieldIndex=21;
SDF_FILE_HDR(21).BinaryIndex=[61 64];
SDF_FILE_HDR(21).FieldName=['offset_of_YDATA_record'];
SDF_FILE_HDR(21).DataType='long';
SDF_FILE_HDR(21).RangeUnits='-1:(2^31)-1';
SDF_FILE_HDR(21).Description = ['Y-axis Data record’s byte offset ' ...
'from beginning of file.'];

SDF_FILE_HDR(22).FieldIndex=22;
SDF_FILE_HDR(22).BinaryIndex=[65 66];
SDF_FILE_HDR(22).FieldName=['num_of_SCAN_BIG_RECORD'];
SDF_FILE_HDR(22).DataType='short';
SDF_FILE_HDR(22).RangeUnits='0:32767';
SDF_FILE_HDR(22).Description = 'total of SDF_SCAN_BIG and SDF_SCAN_VAR records.';

SDF_FILE_HDR(23).FieldIndex=23;
SDF_FILE_HDR(23).BinaryIndex=[67 68];
SDF_FILE_HDR(23).FieldName=['num_of_COMMENT_record'];
SDF_FILE_HDR(23).DataType='short';
SDF_FILE_HDR(23).RangeUnits='0:32767';
SDF_FILE_HDR(23).Description = 'total of SDF_COMMENT_HDR records.';

SDF_FILE_HDR(24).FieldIndex=24;
SDF_FILE_HDR(24).BinaryIndex=[69 72];
SDF_FILE_HDR(24).FieldName=['offset_of_SCAN_BIG_record'];
SDF_FILE_HDR(24).DataType='long';
SDF_FILE_HDR(24).RangeUnits='-1:(2^31)-1';
SDF_FILE_HDR(24).Description = ['the offset (from beginning of file) ' ...
'of the first Scan Big or Scan Variable record.'];

SDF_FILE_HDR(25).FieldIndex=25;
SDF_FILE_HDR(25).BinaryIndex=[73 76];
SDF_FILE_HDR(25).FieldName=['offset_of_next_SDF_FILE'];
SDF_FILE_HDR(25).DataType='long';
SDF_FILE_HDR(25).RangeUnits='-1:(2^31)-1';
SDF_FILE_HDR(25).Description = ['allows more than one logical ' ...
'SDF FILE in a physical file.\n' ...
'This supports multiple independent results taken at the \n' ...
'same time. (For example, a time capture where the span \n' ...
'and center frequencies of each channel are completely \n' ...
'unrelated.) This offset points to the FORMAT_STRUCT \n' ...
'record of the next logical SDF FILE in this physical file. \n' ...
'All offsets in the next SDF FILE are relative to the start of \n' ...
'the FORMAT field of the logical file (that is, "B" \n' ...
'followed by "\\0").'];
end

function  SDF_MEAS_HDR=SDF_MEAS_HDR_Template()
%SDF_MEAS_HDR_TEMPLATE - Provides a empty template containing the requested header information
%
%                 - SDF_MEAS_HDR:    Contains settings of measurement parameters.
%
% Syntax: [retval] = SDF_MEAS_HDR_TEMPLATE()
%
% Inputs: none
%
% Outputs:
%   retval - Extracted SDF header information.
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: N/A
% © Copyright 2016 Commonwealth of Australia, represented by the Department of Defence

SDF_MEAS_HDR(1).FieldIndex=1;
SDF_MEAS_HDR(1).BinaryIndex=[1 2];
SDF_MEAS_HDR(1).FieldName='recordType';
SDF_MEAS_HDR(1).DataType='short';
SDF_MEAS_HDR(1).RangeUnits='11';

SDF_MEAS_HDR(2).FieldIndex=2;
SDF_MEAS_HDR(2).BinaryIndex=[3 6];
SDF_MEAS_HDR(2).FieldName='recordSize';
SDF_MEAS_HDR(2).DataType='long';
SDF_MEAS_HDR(2).RangeUnits='156 bytes';

SDF_MEAS_HDR(3).FieldIndex=3;
SDF_MEAS_HDR(3).BinaryIndex=[7 10];
SDF_MEAS_HDR(3).FieldName='unique_record';
SDF_MEAS_HDR(3).DataType='long';
SDF_MEAS_HDR(3).RangeUnits='-1:0:2^31-1';
SDF_MEAS_HDR(3).Description=['byte offset from the beginning of the file ' ...
'to a record containing an instrument-specific measurement ' ...
'header. This field may be ignored when the file is recalled ' ...
'if it is recalled into an instrument type other than that used ' ...
'to create it.'];

SDF_MEAS_HDR(4).FieldIndex=4;
SDF_MEAS_HDR(4).BinaryIndex=[11 14];
SDF_MEAS_HDR(4).FieldName='centerFreqOld';
SDF_MEAS_HDR(4).DataType='float';
SDF_MEAS_HDR(4).RangeUnits='unit is Hz range is i.d.';
SDF_MEAS_HDR(4).Description=['center frequency. ' ...
'** Prior to version 2.0'];

SDF_MEAS_HDR(5).FieldIndex=5;
SDF_MEAS_HDR(5).BinaryIndex=[15 18];
SDF_MEAS_HDR(5).FieldName='spanFreqOld';
SDF_MEAS_HDR(5).DataType='float';
SDF_MEAS_HDR(5).RangeUnits='unit is Hz range is i.d.';
SDF_MEAS_HDR(5).Description=['frequency span. ' ...
'** Prior to version 2.0'];

SDF_MEAS_HDR(6).FieldIndex=6;
SDF_MEAS_HDR(6).BinaryIndex=[19 22];
SDF_MEAS_HDR(6).FieldName='blockSize';
SDF_MEAS_HDR(6).DataType='long';
SDF_MEAS_HDR(6).RangeUnits='';
SDF_MEAS_HDR(6).Description=['number of time-domain samples taken. This ' ...
'field is only valid for FFT measurements.'];

SDF_MEAS_HDR(7).FieldIndex=7;
SDF_MEAS_HDR(7).BinaryIndex=[23 24];
SDF_MEAS_HDR(7).FieldName='zoomModeOn';
SDF_MEAS_HDR(7).DataType='short';
SDF_MEAS_HDR(7).RangeUnits='0:1';
SDF_MEAS_HDR(7).Description=['zoom mode (0=not zoomed, 1=zoomed). ' ...
'This field is only valid for FFT measurements'];
SDF_MEAS_HDR(7).Table=[...
{'0'} {'non zoomed'};...
{'1'} {'zoomed'};...
];

SDF_MEAS_HDR(8).FieldIndex=8;
SDF_MEAS_HDR(8).BinaryIndex=[25 26];
SDF_MEAS_HDR(8).FieldName='startFreqIndexOld';
SDF_MEAS_HDR(8).DataType='short';
SDF_MEAS_HDR(8).RangeUnits=' 0:last_valid_index(SDF_DATA_HDR)';
SDF_MEAS_HDR(8).Description=['the first alias-protected point on a ' ...
'frequency-domain trace. ' ...
'*** Prior to version 3.0'];

SDF_MEAS_HDR(9).FieldIndex=9;
SDF_MEAS_HDR(9).BinaryIndex=[27 28];
SDF_MEAS_HDR(9).FieldName='stopFreqIndexOld';
SDF_MEAS_HDR(9).DataType='short';
SDF_MEAS_HDR(9).RangeUnits=' 0:last_valid_index(SDF_DATA_HDR)';
SDF_MEAS_HDR(9).Description=['the last alias-protected point on a ' ...
'frequency-domain trace. ' ...
'*** Prior to version 3.0'];

SDF_MEAS_HDR(10).FieldIndex=10;
SDF_MEAS_HDR(10).BinaryIndex=[29 30];
SDF_MEAS_HDR(10).FieldName='averageType';
SDF_MEAS_HDR(10).DataType='short';
SDF_MEAS_HDR(10).RangeUnits='0:6';
SDF_MEAS_HDR(10).Description=[''];
SDF_MEAS_HDR(10).Table=[...
{'0'} {'none'};...
{'1'} {'rms'};...
{'2'} {'rms exponential'};...
{'3'} {'vector'};...
{'4'} {'vector exponential'};...
{'5'} {'continuous peak hold'};...
{'6'} {'peak'};...
];

SDF_MEAS_HDR(11).FieldIndex=11;
SDF_MEAS_HDR(11).BinaryIndex=[31 34];
SDF_MEAS_HDR(11).FieldName='averageNum';
SDF_MEAS_HDR(11).DataType='long';
SDF_MEAS_HDR(11).RangeUnits='range is i.d.';
SDF_MEAS_HDR(11).Description=['number of averages.'];

SDF_MEAS_HDR(12).FieldIndex=12;
SDF_MEAS_HDR(12).BinaryIndex=[35 38];
SDF_MEAS_HDR(12).FieldName='pctOverlap';
SDF_MEAS_HDR(12).DataType='float';
SDF_MEAS_HDR(12).RangeUnits='number between 0 and 1';
SDF_MEAS_HDR(12).Description=['percentage of time-domain samples that are ' ...
'shared between successive time records. This field is only ' ...
'valid for FFT measurements.'];

SDF_MEAS_HDR(13).FieldIndex=13;
SDF_MEAS_HDR(13).BinaryIndex=[39:98];
SDF_MEAS_HDR(13).FieldName='measTitle';
SDF_MEAS_HDR(13).DataType='char[60]';
SDF_MEAS_HDR(13).RangeUnits='i.d.';
SDF_MEAS_HDR(13).Description=['measurement title or automath label.'];

SDF_MEAS_HDR(14).FieldIndex=14;
SDF_MEAS_HDR(14).BinaryIndex=[99 102];
SDF_MEAS_HDR(14).FieldName='videoBandWidth';
SDF_MEAS_HDR(14).DataType='float';
SDF_MEAS_HDR(14).RangeUnits='unit is Hz';
SDF_MEAS_HDR(14).Description=['tells you the bandwidth of the instrument’s ' ...
'video filter. This field is only valid for swept spectrum measurements.'];

SDF_MEAS_HDR(15).FieldIndex=15;
SDF_MEAS_HDR(15).BinaryIndex=[103 110];
SDF_MEAS_HDR(15).FieldName='centerFreq';
SDF_MEAS_HDR(15).DataType='double';
SDF_MEAS_HDR(15).RangeUnits='unit is Hz range is i.d.';
SDF_MEAS_HDR(15).Description=['center frequency'];

SDF_MEAS_HDR(16).FieldIndex=16;
SDF_MEAS_HDR(16).BinaryIndex=[111 118];
SDF_MEAS_HDR(16).FieldName='spanFreq';
SDF_MEAS_HDR(16).DataType='double';
SDF_MEAS_HDR(16).RangeUnits='unit is Hz range is i.d.';
SDF_MEAS_HDR(16).Description=['frequency span'];

SDF_MEAS_HDR(17).FieldIndex=17;
SDF_MEAS_HDR(17).BinaryIndex=[119 126];
SDF_MEAS_HDR(17).FieldName='sweepFreq';
SDF_MEAS_HDR(17).DataType='double';
SDF_MEAS_HDR(17).RangeUnits='unit is Hz range is i.d.';
SDF_MEAS_HDR(17).Description=['current frequency for a swept measurement'];

SDF_MEAS_HDR(18).FieldIndex=18;
SDF_MEAS_HDR(18).BinaryIndex=[127 128];
SDF_MEAS_HDR(18).FieldName='measType';
SDF_MEAS_HDR(18).DataType='short';
SDF_MEAS_HDR(18).RangeUnits='-99:10';
SDF_MEAS_HDR(18).Description=['measurement type'];
SDF_MEAS_HDR(18).Table=[...
{'-99'} {'unknown measurement'};...
{'0'} {'spectrum measurement'};...
{'1'} {'network measurement'};...
{'2'} {'swept measurement'};...
{'3'} {'FFT measurement'};...
{'4'} {'orders measurement'};...
{'5'} {'octave measurement'};...
{'6'} {'capture measurement'};...
{'7'} {'correlation measurement'};...
{'8'} {'histogram measurement'};...
{'9'} {'swept network measurement'};...
{'10'} {'FFT network measurement'};...
];

SDF_MEAS_HDR(19).FieldIndex=19;
SDF_MEAS_HDR(19).BinaryIndex=[129 130];
SDF_MEAS_HDR(19).FieldName='realTime';
SDF_MEAS_HDR(19).DataType='short';
SDF_MEAS_HDR(19).RangeUnits='0:1';
SDF_MEAS_HDR(19).Description=['whether the measurement was continuous in time'];
SDF_MEAS_HDR(19).Table=[...
{'0'} {'non continuous'};...
{'1'} {'continuous'};...
];

SDF_MEAS_HDR(20).FieldIndex=20;
SDF_MEAS_HDR(20).BinaryIndex=[131 132];
SDF_MEAS_HDR(20).FieldName='detection';
SDF_MEAS_HDR(20).DataType='short';
SDF_MEAS_HDR(20).RangeUnits='-99:3';
SDF_MEAS_HDR(20).Description=['detection type'];
SDF_MEAS_HDR(20).Table=[...
{'-99'} {'unknown detection type'};...
{'0'} {'sample detection'};...
{'1'} {'positive peak detection'};...
{'2'} {'negative peak detection'};...
{'3'} {'rose-and-fell detection'};...
];

SDF_MEAS_HDR(21).FieldIndex=21;
SDF_MEAS_HDR(21).BinaryIndex=[133 140];
SDF_MEAS_HDR(21).FieldName='sweepTime';
SDF_MEAS_HDR(21).DataType='double';
SDF_MEAS_HDR(21).RangeUnits='unit is sec range is i.d.';
SDF_MEAS_HDR(21).Description=['actual time for a swept measurement'];

SDF_MEAS_HDR(22).FieldIndex=22;
SDF_MEAS_HDR(22).BinaryIndex=[141 144];
SDF_MEAS_HDR(22).FieldName='startFreqIndex';
SDF_MEAS_HDR(22).DataType='long';
SDF_MEAS_HDR(22).RangeUnits='0:last_valid_index(SDF_DATA_HDR)';
SDF_MEAS_HDR(22).Description=['the first alias-protected point on a ' ...
'frequency-domain trace. This field is only valid for FFT ' ...
'measurements, long'];

SDF_MEAS_HDR(23).FieldIndex=23;
SDF_MEAS_HDR(23).BinaryIndex=[145 148];
SDF_MEAS_HDR(23).FieldName='stopFreqIndex';
SDF_MEAS_HDR(23).DataType='long';
SDF_MEAS_HDR(23).RangeUnits='0:last_valid_index(SDF_DATA_HDR)';
SDF_MEAS_HDR(23).Description=['the last alias-protected point on a ' ...
'frequency-domain trace. This field is only valid for FFT ' ...
'measurements, long'];

SDF_MEAS_HDR(24).FieldIndex=2;
SDF_MEAS_HDR(24).BinaryIndex=[149 156];
SDF_MEAS_HDR(24).FieldName='expAverageNum';
SDF_MEAS_HDR(24).DataType='double';
SDF_MEAS_HDR(24).RangeUnits='range is i.d.';
SDF_MEAS_HDR(24).Description=['—number of exponential averages.'];
end

function  SDF_SCAN_STRUCT=SDF_SCAN_STRUCT_Template()
%SDF_SCAN_STRUCT_HDR_TEMPLATE - Provides a empty template containing the requested header information
%
%                 - SDF_SCAN_STRUCT: Tells you how vectors are organized in the Y-axis Data 
%                                    record when the measurement includes multiple scans of data.
%
% Syntax: [retval] = SDF_SCAN_STRUCT_HDR_TEMPLATE()
%
% Inputs: none
%
% Outputs:
%   retval - Extracted SDF header information.
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: N/A
% © Copyright 2016 Commonwealth of Australia, represented by the Department of Defence

SDF_SCAN_STRUCT(1).FieldIndex=1;
SDF_SCAN_STRUCT(1).BinaryIndex=[1 2];
SDF_SCAN_STRUCT(1).FieldName='recordType';
SDF_SCAN_STRUCT(1).DataType='short';
SDF_SCAN_STRUCT(1).RangeUnits='15';

SDF_SCAN_STRUCT(2).FieldIndex=2;
SDF_SCAN_STRUCT(2).BinaryIndex=[3 6];
SDF_SCAN_STRUCT(2).FieldName='recordSize';
SDF_SCAN_STRUCT(2).DataType='long';
SDF_SCAN_STRUCT(2).RangeUnits='variable';

SDF_SCAN_STRUCT(3).FieldIndex=3;
SDF_SCAN_STRUCT(3).BinaryIndex=[7 8];
SDF_SCAN_STRUCT(3).FieldName='num_of_scan';
SDF_SCAN_STRUCT(3).DataType='short';
SDF_SCAN_STRUCT(3).RangeUnits='1:(215)-1';
SDF_SCAN_STRUCT(3).Description=...
['number of times the instrument collected a complete set of x- and y-axis ' ...
'vectors for all scan-based data types.'];

SDF_SCAN_STRUCT(4).FieldIndex=4;
SDF_SCAN_STRUCT(4).BinaryIndex=[9 10];
SDF_SCAN_STRUCT(4).FieldName='last_scan_index';
SDF_SCAN_STRUCT(4).DataType='short';
SDF_SCAN_STRUCT(4).RangeUnits='0:(num_of_scan)-1';
SDF_SCAN_STRUCT(4).Description=...
['index of the last valid scan'];

SDF_SCAN_STRUCT(5).FieldIndex=5;
SDF_SCAN_STRUCT(5).BinaryIndex=[11 12];
SDF_SCAN_STRUCT(5).FieldName='scan_type';
SDF_SCAN_STRUCT(5).DataType='short';
SDF_SCAN_STRUCT(5).RangeUnits='0:1';
SDF_SCAN_STRUCT(5).Description=...
['tells you how the vectors from different scans ' ...
'are organized in the Y-axis Data record.\n' ...
'0=depth—all scans for the first data type’s vectors ' ...
'followed by all scans for the second data type’s vectors, and so on.\n' ...
'1=scan—all data type’s vectors for the first scan ' ...
'followed by all data type’s vectors for the second scan,and so on.'];
SDF_SCAN_STRUCT(5).Table=[...
{'0'} {'depth'};...
{'1'} {'scan'};...
];

SDF_SCAN_STRUCT(6).FieldIndex=6;
SDF_SCAN_STRUCT(6).BinaryIndex=[13 14];
SDF_SCAN_STRUCT(6).FieldName='scanVar_type';
SDF_SCAN_STRUCT(6).DataType='short';
SDF_SCAN_STRUCT(6).RangeUnits='';
SDF_SCAN_STRUCT(6).Description=...
['tells you the size and format of each scan variable value.\n' ...
'1=short (two-byte, binary-encoded integer)\n' ...
'2=long (four-byte, binary-encoded integer)\n' ...
'3=float (four-byte, binary floating-point number)\n' ...
'4=double (eight-byte, binary floating-point number)'];
SDF_SCAN_STRUCT(6).Table=[...
{'1'} {'short'};...
{'2'} {'long'};...
{'3'} {'float'};...
{'4'} {'double'};...
];

SDF_SCAN_STRUCT(7).FieldIndex=7;
SDF_SCAN_STRUCT(7).BinaryIndex=[15 36];
SDF_SCAN_STRUCT(7).FieldName='scanUnit';
SDF_SCAN_STRUCT(7).DataType='struct';
SDF_SCAN_STRUCT(7).RangeUnits='see SDF_UNIT';
SDF_SCAN_STRUCT(7).Description=...
['engineering unit used for scan variables'];
end

function  SDF_SCANS_BIG =  SDF_SCANS_BIG_Template()
%SDF_SCANS_BIG_TEMPLATE - Provides a empty template containing the requested header information
%
%                 - SDF_SCANS_BIG:   Extended scan header which tells how vectors are organized 
%                                    in the Y-axis data record when the measurement may include 
%                                    more than 32767 scans of data.
%
% Syntax: [retval] = SDF_SCANS_BIG_TEMPLATE()
%
% Inputs: none
%
% Outputs:
%   retval - Extracted SDF header information.
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: N/A
% © Copyright 2016 Commonwealth of Australia, represented by the Department of Defence

SDF_SCANS_BIG(1).FieldIndex=1;
SDF_SCANS_BIG(1).BinaryIndex=[1 2];
SDF_SCANS_BIG(1).FieldName='recordType';
SDF_SCANS_BIG(1).DataType = 'short';
SDF_SCANS_BIG(1).RangeUnits='18';

SDF_SCANS_BIG(2).FieldIndex=2;
SDF_SCANS_BIG(2).BinaryIndex=[3 6];
SDF_SCANS_BIG(2).FieldName='recordSize';
SDF_SCANS_BIG(2).DataType = 'long';
SDF_SCANS_BIG(2).RangeUnits='20';

SDF_SCANS_BIG(3).FieldIndex=3;
SDF_SCANS_BIG(3).BinaryIndex=[7 10];
SDF_SCANS_BIG(3).FieldName='unique_record';
SDF_SCANS_BIG(3).DataType = 'long';
SDF_SCANS_BIG(3).RangeUnits='-1:(231)-1';
SDF_SCANS_BIG(3).Description=...
['byte offset from the beginning of the file to a ' ...
'record containing an instrument-specific scan big header. May ' ...
'be ignored if recalled into a different type instrument.'];

SDF_SCANS_BIG(4).FieldIndex=4;
SDF_SCANS_BIG(4).BinaryIndex=[11 14];
SDF_SCANS_BIG(4).FieldName='num_of_scan';
SDF_SCANS_BIG(4).DataType = 'long';
SDF_SCANS_BIG(4).RangeUnits='-1:(231)-1';
SDF_SCANS_BIG(4).Description=...
['number of times the instrument collected a ' ...
'complete set of x- and y-axis vectors for all scan-based data types.'];

SDF_SCANS_BIG(5).FieldIndex=5;
SDF_SCANS_BIG(5).BinaryIndex=[15 18];
SDF_SCANS_BIG(5).FieldName='last_scan_index';
SDF_SCANS_BIG(5).DataType = 'long';
SDF_SCANS_BIG(5).RangeUnits=':(num_of_scan)-1';
SDF_SCANS_BIG(5).Description=...
['index of the last valid scan.'];

SDF_SCANS_BIG(6).FieldIndex=6;
SDF_SCANS_BIG(6).BinaryIndex=[19 20];
SDF_SCANS_BIG(6).FieldName='scan_type';
SDF_SCANS_BIG(6).DataType = 'short';
SDF_SCANS_BIG(6).RangeUnits='0:1';
SDF_SCANS_BIG(6).Description=...
['tells you how the vectors from different scans ' ...
'are organized in the Y-axis Data record.\n' ...
'0=depth—all scans for the first data type’s vectors ' ...
'followed by all scans for the second data type’s vectors, ' ...
'and so on.\n' ...
'1=scan—all data type’s vectors for the first scan ' ...
'followed by all data type’s vectors for the second scan, ' ...
'and so on.'];
SDF_SCANS_BIG(6).Table=[...
{'0'} {'depth'};...
{'1'} {'scan'};...
];
end

function  SDF_SCANS_VAR=SDF_SCANS_VAR_Template()
%SDF_SCANS_VAR_TEMPLATE - Provides a empty template containing the requested header information
%
%                 - SDF_SCANS_VAR:   Contains a number that identifies every scan 
%                                    (scan time, RPM, or scan number).
%
% Syntax: [retval] = SDF_SCANS_VAR_TEMPLATE()
%
% Inputs: none
%
% Outputs:
%   retval - Extracted SDF header information.
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: N/A
% © Copyright 2016 Commonwealth of Australia, represented by the Department of Defence

SDF_SCANS_VAR(1).FieldIndex=1;
SDF_SCANS_VAR(1).BinaryIndex=[1 2];
SDF_SCANS_VAR(1).FieldName='recordType';
SDF_SCANS_VAR(1).DataType='short';
SDF_SCANS_VAR(1).RangeUnits='19';

SDF_SCANS_VAR(2).FieldIndex=2;
SDF_SCANS_VAR(2).BinaryIndex=[3 6];
SDF_SCANS_VAR(2).FieldName='recordSize';
SDF_SCANS_VAR(2).DataType='long';
SDF_SCANS_VAR(2).RangeUnits='variable';

SDF_SCANS_VAR(3).FieldIndex=3;
SDF_SCANS_VAR(3).BinaryIndex=[7 10];
SDF_SCANS_VAR(3).FieldName='unique_record';
SDF_SCANS_VAR(3).DataType='long';
SDF_SCANS_VAR(3).RangeUnits='-1:(231)-1';
SDF_SCANS_VAR(3).Description=...
['byte offset from the beginning of the file to a ' ...
'record containing an instrument-specific scan variable header. ' ...
'May be ignored if recalled into a different type instrument.'];

SDF_SCANS_VAR(4).FieldIndex=4;
SDF_SCANS_VAR(4).BinaryIndex=[11 14];
SDF_SCANS_VAR(4).FieldName='headersize';
SDF_SCANS_VAR(4).DataType='long';
SDF_SCANS_VAR(4).RangeUnits='54';
SDF_SCANS_VAR(4).Description=...
['size of the header portion of this record ' ...
'(excluding the scan variable values).'];

SDF_SCANS_VAR(5).FieldIndex=5;
SDF_SCANS_VAR(5).BinaryIndex=[15 16];
SDF_SCANS_VAR(5).FieldName='scanBase_type';
SDF_SCANS_VAR(5).DataType='short';
SDF_SCANS_VAR(5).RangeUnits='0:5';
SDF_SCANS_VAR(5).Description=...
['type of scan variable'];
SDF_SCANS_VAR(5).Table=[...
{'0'} {'unknown'};...
{'1'} {'scan number'};...
{'2'} {'time'};...
{'3'} {'RPM'};...
{'4'} {'temperature'};...
{'5'} {'tachometer count'};...
];

SDF_SCANS_VAR(6).FieldIndex=6;
SDF_SCANS_VAR(6).BinaryIndex=[17 18];
SDF_SCANS_VAR(6).FieldName='scanOrder_type';
SDF_SCANS_VAR(6).DataType='short';
SDF_SCANS_VAR(6).RangeUnits='0:2';
SDF_SCANS_VAR(6).Description=...
['progression of scan values'];
SDF_SCANS_VAR(6).Table=[...
{'0'} {'unknown'};...
{'1'} {'increasinhg in value'};...
{'2'} {'decreasing in value'};...
];


SDF_SCANS_VAR(7).FieldIndex=7;
SDF_SCANS_VAR(7).BinaryIndex=[19 20];
SDF_SCANS_VAR(7).FieldName='DATA_recordNum';
SDF_SCANS_VAR(7).DataType='short';
SDF_SCANS_VAR(7).RangeUnits='0:num_of_DATA_HDR_record -1';
SDF_SCANS_VAR(7).Description=...
['SDF_DATA_HDR record number associated with this record, ' ...
'(-1 if no specific association).'];

SDF_SCANS_VAR(8).FieldIndex=8;
SDF_SCANS_VAR(8).BinaryIndex=[21:30];
SDF_SCANS_VAR(8).FieldName='scan_ID';
SDF_SCANS_VAR(8).DataType='char[10]';
SDF_SCANS_VAR(8).RangeUnits='i.d.';
SDF_SCANS_VAR(8).Description=...
['name of scan information.'];

SDF_SCANS_VAR(9).FieldIndex=9;
SDF_SCANS_VAR(9).BinaryIndex=[31 32];
SDF_SCANS_VAR(9).FieldName='scanVar_type';
SDF_SCANS_VAR(9).DataType='short';
SDF_SCANS_VAR(9).RangeUnits='0:4';
SDF_SCANS_VAR(9).Description=...
['tells you the size and format of each scan variable value\n' ...
'1=short (two-byte binary-encoded integer)\n' ...
'2=long (four-byte binary-encoded integer)\n' ...
'3=float (four-byte floating-point number)\n' ...
'4=double (eight-byte floating-point number)'];
SDF_SCANS_VAR(9).Table=[...
{'1'} {'short'};...
{'2'} {'long'};...
{'3'} {'float'};...
{'4'} {'double'};...
];

SDF_SCANS_VAR(10).FieldIndex=10;
SDF_SCANS_VAR(10).BinaryIndex=[33 54];
SDF_SCANS_VAR(10).FieldName='scanUnit';
SDF_SCANS_VAR(10).DataType='struct';
SDF_SCANS_VAR(10).RangeUnits='see SDF_UNIT';
SDF_SCANS_VAR(10).Description=...
['engineering unit used for sca'];
end

function  SDF_UNIT =  SDF_UNIT_Template()
%SDF_UNIT_TEMPLATE - Provides a empty template containing the requested header information
% 
%                 - SDF_UNIT:        Contains eng. units & scaling information for the traces.
%
% Syntax: [retval] = SDF_UNIT_TEMPLATE()
%
% Inputs: none
%
% Outputs:
%   retval - Extracted SDF header information.
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: N/A
% © Copyright 2016 Commonwealth of Australia, represented by the Department of Defence

SDF_UNIT(1).FieldIndex=1;
SDF_UNIT(1).BinaryIndex=[1 10];
SDF_UNIT(1).FieldName='label';
SDF_UNIT(1).DataType = 'char[10]';
SDF_UNIT(1).RangeUnits='i.d.';

SDF_UNIT(2).FieldIndex=2;
SDF_UNIT(2).BinaryIndex=[11 14];
SDF_UNIT(2).FieldName='factor';
SDF_UNIT(2).DataType = 'float';
SDF_UNIT(2).RangeUnits='|10^34(except 0)';


SDF_UNIT(3).FieldIndex=3;
SDF_UNIT(3).BinaryIndex=[15];
SDF_UNIT(3).FieldName='mass';
SDF_UNIT(3).DataType = 'char';
SDF_UNIT(3).RangeUnits='-128:127';
SDF_UNIT(3).Description=...
['2 times the unit exponent for the mass dimension'];

SDF_UNIT(4).FieldIndex=4;
SDF_UNIT(4).BinaryIndex=[16];
SDF_UNIT(4).FieldName='length';
SDF_UNIT(4).DataType = 'char';
SDF_UNIT(4).RangeUnits='-128:127';
SDF_UNIT(4).Description=...
['2 times the unit exponent for the length dimension'];

SDF_UNIT(5).FieldIndex=5;
SDF_UNIT(5).BinaryIndex=[17];
SDF_UNIT(5).FieldName='time';
SDF_UNIT(5).DataType = 'char';
SDF_UNIT(5).RangeUnits='-128:127';
SDF_UNIT(5).Description=...
['2 times the unit exponent for the time dimension'];

SDF_UNIT(6).FieldIndex=6;
SDF_UNIT(6).BinaryIndex=[18];
SDF_UNIT(6).FieldName='current';
SDF_UNIT(6).DataType = 'char';
SDF_UNIT(6).RangeUnits='-128:127';
SDF_UNIT(6).Description=...
['2 times the unit exponent for the current dimension'];

SDF_UNIT(7).FieldIndex=7;
SDF_UNIT(7).BinaryIndex=[19];
SDF_UNIT(7).FieldName='temperature';
SDF_UNIT(7).DataType = 'char';
SDF_UNIT(7).RangeUnits='-128:127';
SDF_UNIT(7).Description=...
['2 times the unit exponent for the temperature dimension'];

SDF_UNIT(8).FieldIndex=8;
SDF_UNIT(8).BinaryIndex=[20];
SDF_UNIT(8).FieldName='luminal_intensity';
SDF_UNIT(8).DataType = 'char';
SDF_UNIT(8).RangeUnits='-128:127';
SDF_UNIT(8).Description=...
['2 times the unit exponent for the luminal intensity dimension'];

SDF_UNIT(9).FieldIndex=9;
SDF_UNIT(9).BinaryIndex=[21];
SDF_UNIT(9).FieldName='mole';
SDF_UNIT(9).DataType = 'char';
SDF_UNIT(9).RangeUnits='-128:127';
SDF_UNIT(9).Description=...
['2 times the unit exponent for the mole dimension'];

SDF_UNIT(10).FieldIndex=10;
SDF_UNIT(10).BinaryIndex=[22];
SDF_UNIT(10).FieldName='plane_angle';
SDF_UNIT(10).DataType = 'char';
SDF_UNIT(10).RangeUnits='-128:127';
SDF_UNIT(10).Description=...
['2 times the unit exponent for the plane angle dimension'];
SDF_UNIT(10).Table={};
end

function  SDF_VECTOR_HDR=SDF_VECTOR_HDR_Template()
%SDF_VECTOR_HDR_TEMPLATE - Provides a empty template containing the requested header information
%
%                 - SDF_VECTOR_HDR:  Tells you which channel (or pair of channels) provided data 
%                                    for a single trace.
%
% Syntax: [retval] = SDF_VECTOR_HDR_TEMPLATE()
%
% Inputs: none
%
% Outputs:
%   retval - Extracted SDF header information.
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: N/A
% © Copyright 2016 Commonwealth of Australia, represented by the Department of Defence

SDF_VECTOR_HDR(1).FieldIndex=1;
SDF_VECTOR_HDR(1).BinaryIndex=[1 2];
SDF_VECTOR_HDR(1).FieldName='recordType';
SDF_VECTOR_HDR(1).DataType='short';
SDF_VECTOR_HDR(1).RangeUnits='13';

SDF_VECTOR_HDR(2).FieldIndex=2;
SDF_VECTOR_HDR(2).BinaryIndex=[3 6];
SDF_VECTOR_HDR(2).FieldName='recordSize';
SDF_VECTOR_HDR(2).DataType='long';
SDF_VECTOR_HDR(2).RangeUnits='18 bytes';

SDF_VECTOR_HDR(3).FieldIndex=3;
SDF_VECTOR_HDR(3).BinaryIndex=[7 10];
SDF_VECTOR_HDR(3).FieldName='unique_record';
SDF_VECTOR_HDR(3).DataType='long';
SDF_VECTOR_HDR(3).RangeUnits='-1:0:2^31-1';
SDF_VECTOR_HDR(3).Description=...
['byte offset from the beginning of the file to ' ...
'a record containing an instrument-specific vector header. ' ...
'May be ignored if recalled into a different type instrument.'];

SDF_VECTOR_HDR(4).FieldIndex=4;
SDF_VECTOR_HDR(4).BinaryIndex=[11 12];
SDF_VECTOR_HDR(4).FieldName='the_CHANNEL_record(1)';
SDF_VECTOR_HDR(4).DataType='short';
SDF_VECTOR_HDR(4).RangeUnits='-1, 0:32767';
SDF_VECTOR_HDR(4).Description=...
['tells you which channel or channels provided data for this trace.\n' ...
'Each element of this array contains an index to a Channel Header record.\n' ...
'An element refers to no channel if the index value is -1.'];

SDF_VECTOR_HDR(6).FieldIndex=4;
SDF_VECTOR_HDR(6).BinaryIndex=[13 14];
SDF_VECTOR_HDR(6).FieldName='the_CHANNEL_record(2)';
SDF_VECTOR_HDR(6).DataType='short';
SDF_VECTOR_HDR(6).RangeUnits='-1, 0:32767';
SDF_VECTOR_HDR(6).Description=...
['tells you which channel or channels provided data for this trace.\n' ...
'Each element of this array contains an index to a Channel Header record.\n' ...
'An element refers to no channel if the index value is -1.'];

SDF_VECTOR_HDR(5).FieldIndex=5;
SDF_VECTOR_HDR(5).BinaryIndex=[15 16];
SDF_VECTOR_HDR(5).FieldName='powrOfChan(1)';
SDF_VECTOR_HDR(5).DataType='short';
SDF_VECTOR_HDR(5).RangeUnits='0:32767';
SDF_VECTOR_HDR(5).Description=...
['tells you what exponent was applied to corresponding channel data to ' ...
'create this trace. (For example, pwrOfChan[0] is the exponent applied to ' ...
'the_CHANNEL record[0]’s data.)\n\n' ...
'record [0] for row channel\n' ...
'Chan [0] “ ��? \n\n' ...
'record [1] for column channel\n' ...
'Chan [1] “ ��? \n\n' ...
'NOTE: You must divide each pwrOfChan value by 48 to obtain ' ...
'the true value of the exponent.'];

SDF_VECTOR_HDR(7).FieldIndex=5;
SDF_VECTOR_HDR(7).BinaryIndex=[17 18];
SDF_VECTOR_HDR(7).FieldName='powrOfChan(2)';
SDF_VECTOR_HDR(7).DataType='short';
SDF_VECTOR_HDR(7).RangeUnits='0:32767';
SDF_VECTOR_HDR(7).Description=...
['tells you what exponent was applied to corresponding channel data to ' ...
'create this trace. (For example, pwrOfChan[0] is the exponent applied to ' ...
'the_CHANNEL record[0]’s data.)\n\n' ...
'record [0] for row channel\n' ...
'Chan [0] “ ��? \n\n' ...
'record [1] for column channel\n' ...
'Chan [1] “ ��? \n\n' ...
'NOTE: You must divide each pwrOfChan value by 48 to obtain ' ...
'the true value of the exponent.'];
SDF_VECTOR_HDR(7).Table={};
end

function  SDF_WINDOW =  SDF_WINDOW_Template()
%SDF_WINDOW_TEMPLATE - Provides a empty template containing the requested header information
% 
%                 - SDF_WINDOW:      Contains windowing information for frequency domain traces.
%
% Syntax: [retval] = SDF_WINDOW_TEMPLATE()
%
% Inputs: none
%
% Outputs:
%   retval - Extracted SDF header information.
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: N/A
% © Copyright 2016 Commonwealth of Australia, represented by the Department of Defence

SDF_WINDOW(1).FieldIndex=1;
SDF_WINDOW(1).BinaryIndex=[1 2];
SDF_WINDOW(1).FieldName='windowType';
SDF_WINDOW(1).DataType = 'short';
SDF_WINDOW(1).RangeUnits='0:12';
SDF_WINDOW(1).Table=[...
{'0'} {'window not applied'};...
{'1'} {'Hanning'};...
{'2'} {'Flat Top'};...
{'3'} {'Uniform'};...
{'4'} {'Force'};...
{'5'} {'Response'};...
{'6'} {'user-defined'};...
{'7'} {'Hamming'};...
{'8'} {'P301'};...
{'9'} {'P310'};...
{'10'} {'Kaiser-Bessel'};...
{'11'} {'Harris'};...
{'12'} {'Blackman'};...
{'13'} {'Resolution filter'};...
{'14'} {'Correlation Lead Lag'};...
{'15'} {'Correlation Lag'};...
{'16'} {'Gated'};...
{'17'} {'P400'};...
]; 

SDF_WINDOW(2).FieldIndex=2;
SDF_WINDOW(2).BinaryIndex=[3 4];
SDF_WINDOW(2).FieldName='windowCorrMode';
SDF_WINDOW(2).DataType = 'short';
SDF_WINDOW(2).RangeUnits='0:2';
SDF_WINDOW(2).Table=[...
{'0'} {'correction not applied'};...
{'1'} {'narrow band correction applied'};...
{'2'} {'wide band correction applied'};...
];

SDF_WINDOW(3).FieldIndex=3;
SDF_WINDOW(3).BinaryIndex=[5 8];
SDF_WINDOW(3).FieldName='windowBandWidth';
SDF_WINDOW(3).DataType = 'float';
SDF_WINDOW(3).RangeUnits='unit is bins unit is Hz (if windowType=13), range is i.d.';
SDF_WINDOW(3).Descripton=...
['NOTE: When windowType = 13, this field contains the ' ...
'instrument’s resolution bandwidth.'];

SDF_WINDOW(4).FieldIndex=4;
SDF_WINDOW(4).BinaryIndex=[9 12];
SDF_WINDOW(4).FieldName='windowTimeConst';
SDF_WINDOW(4).DataType = 'long';
SDF_WINDOW(4).RangeUnits='unit is sec, range is i.d.';
SDF_WINDOW(4).Descripton=...
['determines decay of Force and Response windows'];

SDF_WINDOW(5).FieldIndex=5;
SDF_WINDOW(5).BinaryIndex=[13 16];
SDF_WINDOW(5).FieldName='windowTrunc';
SDF_WINDOW(5).DataType = 'float';
SDF_WINDOW(5).RangeUnits='unit is sec, range is i.d.';
SDF_WINDOW(5).Descripton=...
['width of FORCE window'];

SDF_WINDOW(6).FieldIndex=6;
SDF_WINDOW(6).BinaryIndex=[17 20];
SDF_WINDOW(6).FieldName='wideBandCorr';
SDF_WINDOW(6).DataType = 'float';
SDF_WINDOW(6).RangeUnits='|10^34 (except 0)';
SDF_WINDOW(6).Descripton=...
['correction factor for wide-band signals (like random noise)'];

SDF_WINDOW(7).FieldIndex=7;
SDF_WINDOW(7).BinaryIndex=[21 24];
SDF_WINDOW(7).FieldName='narrowBandCorr';
SDF_WINDOW(7).DataType = 'float';
SDF_WINDOW(7).RangeUnits='|10^34 (except 0)';
SDF_WINDOW(7).Descripton=...
['correction factor for narrow band signals (like sinesoidal wave)'];
end

function  SDF_XDATA_HDR =  SDF_XDATA_HDR_Template()
%SDF_XDATA_HDR_TEMPLATE - Provides a empty template containing the requested header information
%
%                 - SDF_XDATA_HDR:   Contains the x-axis data needed to reconstruct any trace.
%
% Syntax: [retval] = SDF_XDATA_HDR_TEMPLATE()
%
% Inputs: none
%
% Outputs:
%   retval - Extracted SDF header information.
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: N/A
% © Copyright 2016 Commonwealth of Australia, represented by the Department of Defence

SDF_XDATA_HDR(1).FieldIndex=1;
SDF_XDATA_HDR(1).BinaryIndex=[1 2];
SDF_XDATA_HDR(1).FieldName='recordType';
SDF_XDATA_HDR(1).DataType = 'short';
SDF_XDATA_HDR(1).RangeUnits='20';

SDF_XDATA_HDR(2).FieldIndex=2;
SDF_XDATA_HDR(2).BinaryIndex=[3 6];
SDF_XDATA_HDR(2).FieldName='recordSize';
SDF_XDATA_HDR(2).DataType = 'long';
SDF_XDATA_HDR(2).RangeUnits='variable';
SDF_XDATA_HDR(2).Table = {};
end

function  SDF_YDATA_HDR =  SDF_YDATA_HDR_Template()
%SDF_YDATA_HDR_TEMPLATE - Provides a empty template containing the requested header information
%
%                 - SDF_YDATA_HDR:   Contains the y-axis data needed to reconstruct any trace.
%
% Syntax: [retval] = SDF_YDATA_HDR_TEMPLATE()
%
% Inputs: none
%
% Outputs:
%   retval - Extracted SDF header information.
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: N/A
% © Copyright 2016 Commonwealth of Australia, represented by the Department of Defence

SDF_YDATA_HDR(1).FieldIndex=1;
SDF_YDATA_HDR(1).BinaryIndex=[1 2];
SDF_YDATA_HDR(1).FieldName='recordType';
SDF_YDATA_HDR(1).DataType = 'short';
SDF_YDATA_HDR(1).RangeUnits='20';

SDF_YDATA_HDR(2).FieldIndex=2;
SDF_YDATA_HDR(2).BinaryIndex=[3 6];
SDF_YDATA_HDR(2).FieldName='recordSize';
SDF_YDATA_HDR(2).DataType = 'long';
SDF_YDATA_HDR(2).RangeUnits='variable';
SDF_YDATA_HDR(2).Table = {};
end

function  SDF_UNIQUE =  SDF_UNIQUE_Template()
%SDF_UNIQUE_HDR_TEMPLATE - Provides a empty template containing the requested header information
% 
%                 - SDF_UNIQUE:      Makes the SDF file format flexible. The eight common record
%                                    types define parameters that are common to many instruments 
%                                    and systems. However, a particular instrument or system may 
%                                    need to save and recall the states of additional parameters. 
%                                    These states can reside in Unique records.
%
% Syntax: [retval] = SDF_UNIQUE_HDR_TEMPLATE()
%
% Inputs: none
%
% Outputs:
%   retval - Extracted SDF header information.
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Justin Dinale
% DST Group, Department of Defence
% email: Justin.Dinale@dst.defence.gov.au
% Website: http://www.dst.defence.gov.au
% November 2016; Latest Version: N/A
% © Copyright 2016 Commonwealth of Australia, represented by the Department of Defence

SDF_UNIQUE(1).FieldIndex=1;
SDF_UNIQUE(1).BinaryIndex=[1 2];
SDF_UNIQUE(1).FieldName='recordType';
SDF_UNIQUE(1).DataType = 'short';
SDF_UNIQUE(1).RangeUnits='20';

SDF_UNIQUE(2).FieldIndex=2;
SDF_UNIQUE(2).BinaryIndex=[3 6];
SDF_UNIQUE(2).FieldName='recordSize';
SDF_UNIQUE(2).DataType = 'long';
SDF_UNIQUE(2).RangeUnits='variable';
SDF_UNIQUE(2).Table = {};
end

