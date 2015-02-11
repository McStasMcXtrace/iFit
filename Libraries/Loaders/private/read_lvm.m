function data = read_lvm(filename,verbose)
%LVM_IMPORT Imports data from a LabView LVM file
% DATA = LVM_IMPORT(FILENAME,VERBOSE) returns the data from a LVM (.lvm)
%  ASCII text file created by LabView.
%
% FILENAME    The name of the .lvm file, with or without ".lvm" extension
%
% VERBOSE     How many messages to display. Default is 1 (few messages),
%              0 = silent, 2 = display file information and all messages
%
% DATA        The data found in the LVM file. DATA is a structure with 
%              fields corresponding to the Segments in the file (see below) 
%              and LVM file header information.
%
%
% This function imports data from a text-formatted LabView Measurement File
%  (LVM, extension ".lvm") into MATLAB. A LVM file can have multiple
%  Segments, so that multiple measurements can be combined in a single
%  file. The output variable DATA is a structure with fields named
%  'Segment1', 'Segment2', etc. Each Segment field is a structure with
%  details about the data in the Segment and the actual data in the field
%  named 'data'. The column labels and units are stored as cell arrays that
%  correspond to the columns in the array of data.
% The size of the data array depends on the type of x-axis data that is
%  stored in the LVM file and the number of channels (num_channels).
%  There are three cases:
%  1) No x-data is included in the file ('No')
%   The data array will have num_channels columns (one column per channel
%   of data).
%  2) One column of x-data is included in the file ('One')
%   The first column of the data array will be the x-values, and the data
%   array will have num_channels+1 columns.
%  3) Each channel has its own x-data ('Multi')
%   Each channel has two columns, one for x-values, and one for data. The
%   data array will have num_channels*2 columns, with the x-values and
%   corresponding data in alternating columns. For example, in a Segment
%   with 4 channels, columns 1,3,5,7 will be the x-values for the data in
%   columns 2,4,6,8.
%
% Note: because MATLAB only works with a "." decimal separator, importing
%  large LVM files that use a "," (or other character) will be noticeably
%  slower. Use a "." decimal separator to avoid this issue.
%
% The LVM file specification is available at:
%   http://zone.ni.com/devzone/cda/tut/p/id/4139
%
% Contrib from:
% http://www.mathworks.com/matlabcentral/fileexchange/19913-lvm-file-import
%
%
% Example:
%
%  Use the following command to read in the data from a file containing two
%   Segments:
%
% >> d=lvm_import('testfile.lvm');
%
% Importing testfile.lvm:
%
% Import complete. 2 Segments found.
%
% >> d
% d = 
%       X_Columns: 'One'
%            user: 'hopcroft'
%     Description: 'Pressure, Flowrate, Heat, Power, Analog Voltage, Pump on, Temp'
%            date: '2008/03/26'
%            time: '12:18:02.156616'
%           clock: [2008 3 26 12 18 2.156616]
%        Segment1: [1x1 struct]
%        Segment2: [1x1 struct]
%
% >> d.Segment1
% ans = 
%            Notes: 'Some notes regarding this data set'
%     num_channels: 8
%          y_units: {8x1 cell}
%          x_units: {8x1 cell}
%               X0: [8x1 double]
%          Delta_X: [8x1 double]
%    column_labels: {9x1 cell}
%             data: [211x9 double]
%          Comment: 'This data rulz'
%
% >> d.Segment1.column_labels{2}
% ans =
% Thermocouple1
%
% >> plot(d.Segment1.data(:,1),d.Segment1.data(:,2));
% >> xlabel(d.Segment1.column_labels{1});
% >> ylabel(d.Segment1.column_labels{2});
%
%
%
% M.A. Hopcroft
%      < mhopeng at gmail.com >
%

% MH Sep2013
% v2.2  fixes for case of comma separator in multi-segment files
%       use cell2mat for performance improvement
%       (thanks to <die-kenny@t-online.de> for bug report and testing)
% MH May2012
% v2.1  handle "no separator" bug
%       (thanks to <adnan.cheema@gmail.com> for bug report and testing)
%       code & comments cleanup
%       remove extraneous column labels (X_Value for "No X" files; Comment)
%       clean up verbose output
%       change some field names to NI names ("Delta_X","X_Columns","Date")
% MH Mar2012
% v2.0  fix "string bug" related to comma-separated decimals 
%       handle multiple Special Headers correctly
%       fix help comments
%       increment version number to match LabView LVM writer
% MH Sep2011
% v1.3  handles LVM Writer version 2.0 (files with decimal separator)
%       Note: if you want to work with older files with a non-"." decimal
%       separator character, change the value of "data.Decimal_Separator"
% MH Sep2010
% v1.2  bugfixes for "Special" header in LVM files.
%        (Thanks to <bobbyjoe23928@gmail.com> for suggestions)
% MH Apr2010
% v1.1  use case-insensitive comparisons to maintain compatibility with
%        NI LVM Writer version 1.00
%
% MH MAY2009
% v1.02 Add filename input
% MH SEP2008
% v1.01 Fix comments, add Cells
% v1.00 Handle all three possibilities for X-columns (No,One,Multi)
%       Handle LVM files with no header
% MH AUG2008
% v0.92 extracts Comment for each Segment
% MH APR2008
% v0.9  initial version
%

%  M. A. Hopcroft, 13 May 2008 (Updated 23 Sep 2013) 
% <http://www.mathworks.com/matlabcentral/fileexchange/19913-lvm-file-import>
%
%Copyright (c) 2013, M. A. Hopcroft
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are
%met:
%
%* Redistributions of source code must retain the above copyright
%  notice, this list of conditions and the following disclaimer.
%* Redistributions in binary form must reproduce the above copyright
%  notice, this list of conditions and the following disclaimer in
%  the documentation and/or other materials provided with the distribution
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

%#ok<*ASGLU>

% message level
if nargin < 2, verbose = 0; end
if verbose >= 1, fprintf(1,'\nlvm_import v2.2\n'); end

% ask for filename if not provided already
if nargin < 1
    filename=input(' Enter the name of the .lvm file: ','s');
    fprintf(1,'\n');
end


%% Open the data file
% open and verify the file
fid=fopen(filename);
if fid ~= -1, % then file exists
    fclose(fid);
else
    filename=strcat(filename,'.lvm');
    fid=fopen(filename);
    if fid ~= -1, % then file exists
        fclose(fid);
    else
        error(['File not found in current directory! (' pwd ')']);
    end
end

% open the validated file
fid=fopen(filename);

if verbose >= 1, fprintf(1,' Importing %s:\n\n',filename); end
if verbose >= 2, fprintf(1,' File Header:\n'); end

% is it really a LVM file?
linein=fgetl(fid);
if verbose >= 2, fprintf(1,'%s\n',linein); end
if ~strcmp(sscanf(linein,'%s'),'LabVIEWMeasurement')
    try
        data.Segment1.data = dlmread(filename,'\t');
        if verbose >= 1, fprintf(1,'This file appears to be an LVM file with no header.\n'); end
        if verbose >= 1, fprintf(1,'Data was copied, but no other information is available.\n'); end
        return
    catch fileEx
        error([ mfilename ': This does not appear to be a text-format LVM file (no header).' ]);
    end
end


%% Process file header
% The file header contains several fields with useful information

% default values
data.Decimal_Separator = '.';
text_delimiter='\t';
data.X_Columns='One';

% File header contains date, time, etc.
% Also the file delimiter and decimal separator (LVM v2.0)
while 1 
    
    % get a line from the file
    linein=fgetl(fid);
    % handle spurious carriage returns
    if isempty(linein), linein=fgetl(fid); end
    if verbose >= 2, fprintf(1,'%s\n',linein); end
    % what is the tag for this line?
    t_in = textscan(linein,'%s');
    if isempty(t_in{1})
        tag='notag';
    else
        tag = t_in{1}{1};
    end
    % exit when we reach the end of the header
    if strcmpi(tag,'***End_of_Header***')
        if verbose >= 2, fprintf(1,'\n'); end
        break
    end
    
    % get the value corresponding to the tag
    if ~strcmp(tag,'notag')
        v_in = textscan(linein,'%*s %s','delimiter','\t','whitespace','','MultipleDelimsAsOne', 1);
        if ~isempty(v_in{1})
            val = v_in{1}{1};
    
            switch tag
                case 'Date'
                    data.Date = val;
                case 'Time'
                    data.Time = val;
                case 'Operator'
                    data.user = val;
                case 'Description'
                    data.Description = val;
                case 'Project'
                    data.Project = val;            
                case 'Separator'
                    if strcmp(val,'Tab')
                        text_delimiter='\t';
                    elseif strcmp(val,'Comma')
                        text_delimiter=',';
                    end
                case 'X_Columns'
                    data.X_Columns = val;
                case 'Decimal_Separator'
                    data.Decimal_Separator = val;
            end
            
        end
    end    
    
end

% create matlab-formatted date vector
if isfield(data,'time') && isfield(data,'date')
    dt = textscan(data.Date,'%d','delimiter','/');
    tm = textscan(data.Time,'%d','delimiter',':');
    if length(tm{1})==3
        data.clock=[dt{1}(1) dt{1}(2) dt{1}(3) tm{1}(1) tm{1}(2) tm{1}(3)];
    elseif length(tm{1})==2
        data.clock=[dt{1}(1) dt{1}(2) dt{1}(3) tm{1}(1) tm{1}(2) 0];
    else
        data.clock=[dt{1}(1) dt{1}(2) dt{1}(3) 0 0 0];
    end
end



%% Process segments
% process data segments in a loop until finished
segnum = 1;

while 1
    %segnum = segnum +1;
    fieldnm = ['Segment' num2str(segnum)];

    %% - Segment header
    if verbose >= 1, fprintf(1,' Segment %d:\n\n',segnum); end
    % loop to read segment header
    while 1
        % get a line from the file
        linein=fgetl(fid);
        % handle spurious carriage returns/blank lines/end of file
        while isempty(linein), linein=fgetl(fid); end
        if feof(fid), break; end
        if verbose >= 2, fprintf(1,'%s\n',linein); end
        
        % Ignore "special segments"
        % "special segments" can hold other types of data. The type tag is
        % the first line after the Start tag. As of version 2.0,
        % LabView defines three types:
        %  Binary_Data
        %  Packet_Notes
        %  Wfm_Sclr_Meas
        % In theory, users can define their own types as well. LVM_IMPORT
        %  ignores any "special segments" it finds.
        % If special segments are handled in future versions, recommend
        %  moving the handler outside the segment read loop.
        if strfind(linein,'***Start_Special***')            
            special_seg = 1;
            while special_seg                

                while 1 % process lines until we find the end of the special segment 
                    % get a line from the file
                    linein=fgetl(fid);
                    % handle spurious carriage returns
                    if isempty(linein), linein=fgetl(fid); end
                    % test for end of file
                    if linein==-1, break; end
                    if verbose >= 2, fprintf(1,'%s\n',linein); end
                    if strfind(linein,'***End_Special***')
                        if verbose >= 2, fprintf(1,'\n'); end
                        break
                    end
                end
                
                % get the next line and proceed with file 
                %  (there may be additional Special Segments)
                linein=fgetl(fid);
                % handle spurious carriage returns/blank lines/end of file
                while isempty(linein), linein=fgetl(fid); end
                if feof(fid), break; end
                if isempty(strfind(linein,'***Start_Special***'))
                    special_seg = 0;
                    if verbose >= 1, fprintf(1,' [Special Segment ignored]\n\n'); end
                end
            end
        end % end special segment handler
        
        
        % what is the tag for this line?
        t_in = textscan(linein,'%s');
        if isempty(t_in{1})
            tag='notag';
        else
            tag = t_in{1}{1};
        end
        % exit when we reach the end of the header
        if strcmpi(tag,'***End_of_Header***')
            if verbose >= 2, fprintf(1,'\n'); end
            break
        end
   
        
        switch tag
            case 'Notes'
                %d_in = textscan(linein,'%*s %s','delimiter','\t','whitespace','');
                d_in = linein;
                data.(fieldnm).Notes=d_in;
            case 'Test_Name'
                %d_in = textscan(linein,'%*s %s','delimiter','\t','whitespace','');
                d_in = linein;
                data.(fieldnm).Test_Name = d_in;  %d_in{1}{1};           
            case 'Channels'
                numchan = textscan(linein,'%*s %d',1);
                data.(fieldnm).num_channels = numchan{1};
            case 'Samples'
                numsamp = textscan(linein,'%s','delimiter',text_delimiter);
                numsamp1 = numsamp{1};
                numsamp1(1)=[]; % remove tag "Samples"
                numsamp2=str2num(cell2mat(numsamp1));                              %#ok<ST2NM>
                data.(fieldnm).num_samples = numsamp2(:)';                
            case 'Y_Unit_Label'
                Y_units = textscan(linein,'%s','delimiter',text_delimiter);
                data.(fieldnm).y_units=Y_units{1}';
                data.(fieldnm).y_units(1)=[]; % remove tag
            case 'Y_Dimension'
                Y_Dim = textscan(linein,'%s','delimiter',text_delimiter);
                data.(fieldnm).y_type=Y_Dim{1}';
                data.(fieldnm).y_type(1)=[]; % remove tag
            case 'X_Unit_Label'
                X_units = textscan(linein,'%s','delimiter',text_delimiter);
                data.(fieldnm).x_units=X_units{1}';
                data.(fieldnm).x_units(1)=[];
            case 'X_Dimension'
                X_Dim = textscan(linein,'%s','delimiter',text_delimiter);
                data.(fieldnm).x_type=X_Dim{1}';
                data.(fieldnm).x_type(1)=[]; % remove tag
            case 'X0'           
                [Xnought, val]=strtok(linein);
                if ~strcmp(data.Decimal_Separator,'.')
                    val = strrep(val,data.Decimal_Separator,'.');
                end
                data.(fieldnm).X0 = sscanf(val,'%e');
            case 'Delta_X' %,
                [Delta_X, val]=strtok(linein);
                if ~strcmp(data.Decimal_Separator,'.')
                    val = strrep(val,data.Decimal_Separator,'.');
                end
                data.(fieldnm).Delta_X = sscanf(val,'%e');                
        end
        
    end % end reading segment header loop
    

    % after each segment header is the row of column labels
    linein=fgetl(fid);
    Y_labels = textscan(linein,'%s','delimiter',text_delimiter);       
    data.(fieldnm).column_labels=Y_labels{1}';
    % The X-column always exists, even if it is empty. Remove if not used.
    if strcmpi(data.X_Columns,'No')
        data.(fieldnm).column_labels(1)=[];
    end
    % remove empty entries and "Comment" label
    if any(strcmpi(data.(fieldnm).column_labels,'Comment'))
        data.(fieldnm).column_labels=data.(fieldnm).column_labels(1:find(strcmpi(data.(fieldnm).column_labels,'Comment'))-1);
    end
    % display column labels
    if verbose >= 1
        fprintf(1,' Data Columns:\n | ');
        for i=1:length(data.(fieldnm).column_labels)
            fprintf(1,'%s | ',data.(fieldnm).column_labels{i});
        end
        fprintf(1,'\n\n');
    end



    %% - Segment Data
    % Create a format string for textscan depending on the number/type of
    %  channels. If there are additional segments, texscan will quit when
    %  it comes to a text line which does not fit the format, and the loop
    %  will repeat.
    if verbose >= 1, fprintf(1,' Importing data from Segment %d...',segnum); end
    
    % How many data columns do we have? (including X data)
    switch data.X_Columns
        case 'No'
            % an empty X column exists in the file
            numdatacols = data.(fieldnm).num_channels+1;
            xColPlural='no X-Columns';
        case 'One'
            numdatacols = data.(fieldnm).num_channels+1;
            xColPlural='one X-Column';
        case 'Multi'
            numdatacols = data.(fieldnm).num_channels*2;
            xColPlural='multiple X-Columns';
    end
    
    % handle case of not using periods (aka "dot" or ".") for decimal point separators
    %  (LVM version 2.0+)
    if ~strcmp(data.Decimal_Separator,'.')
        if verbose >= 2, fprintf(1,'\n  (using decimal separator "%s")\n',data.Decimal_Separator); end
        
        % create a format string for reading data as numbers
        fs = '%f'; for i=2:numdatacols, fs = [fs ' %f']; end                %#ok<AGROW>
        % add one more column for the comment field
        fs = [fs ' %s'];                                                   %#ok<AGROW>
        rawdata=[];
        while 1
            cline=fgetl(fid);
            if isempty(cline) || (isnumeric(cline) && cline==-1), break; end
            cline=strrep(cline,',','.');
            cread = textscan(cline,fs,'delimiter',text_delimiter);
            if isempty(rawdata) && all(size(cread{size(cread,2)})) > 0
                % save first row comment as The Comment for this segment
                data.(fieldnm).Comment = cread{size(cread,2)}{1};
            end
            rawdata=[rawdata; cread(1:numdatacols)];                   %#ok<AGROW>
        end

    else
        % create a format string for reading data as numbers
        fs = '%f'; for i=2:numdatacols, fs = [fs ' %f']; end                    %#ok<AGROW>
        % add one more column for the comment field
        fs = [fs ' %s'];                                                        %#ok<AGROW>
        % read the data from file
        rawdata = textscan(fid,fs,'delimiter',text_delimiter);
        % save first row comment as The Comment for this segment
        data.(fieldnm).Comment = rawdata{size(rawdata,2)}{1};
    end
    
    % v2.2 use cell2mat here instead of a loop for better performance
    % consolidate data into a simple array, ignore comments
    data.(fieldnm).data=cell2mat(rawdata(:,1:numdatacols));
    
    % If we have a "No X data" file, remove the first column (it is empty/NaN)
    if strcmpi(data.X_Columns,'No')
        data.(fieldnm).data=data.(fieldnm).data(:,2:end);
    end
    
    if verbose >= 1, fprintf(1,' complete (%g data points).\n\n',length(data.(fieldnm).data)); end
    
    % test for end of file
    if feof(fid)
        if verbose >= 2, fprintf(1,' [End of File]\n\n'); end
        break;
    else
        segnum = segnum+1;
    end    
    
    
end % end process segment


if verbose >= 1, fprintf(1,'\n Import complete. File has %s and %d Data Segments.\n\n',xColPlural,segnum); end


% close the file
fclose(fid);
return
