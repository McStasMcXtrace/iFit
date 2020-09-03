function allData = tsread(filename)
% TSREAD reads in parameters from suitable Bruker NMR files
% output = tsread(filename) reads in the contents of the JCAMP-format
% Bruker TopSpin parameter files e.g. acqu, proc etc., and then handles the
% custom fields so they give sensible values. Returns a structure.

% <https://fr.mathworks.com/matlabcentral/fileexchange/71731-tsread>
%  Copyright (c) 2019, Geoffrey Akien
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are met:
%
%  * Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
%
%  * Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution
%  * Neither the name of Lancaster University nor the names of its
%    contributors may be used to endorse or promote products derived from this
%    software without specific prior written permission.
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
%  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
%  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
%  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
%  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
%  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% checks the number of arguments
if nargin ~=1, error([ mfilename ': require 1 input arg' ]); end

% read in all the data as parameter-value pairs
allData = jcampread2(filename);
parameters = allData.Notes(:, 1);
values = allData.Notes(:, 2);

% preparation for the loop
expectedArraySize = zeros(size(parameters));
output = cell(size(parameters, 1), 3);
counter = 0;

% cleanup for multi-line fields first
for m = 1:size(parameters, 1)
    
    if strcmp(parameters{m}, 'END')
        % exception - we don't need to record this, and this has an empty
        % value anyway
        continue;
        
    elseif isempty(values{m})
        % then its part of a multi-line field (where the data is actually
        % stored in parameters, and we have to assume its part of the
        % values from the previous one
        
        % need to parse the previous field
        % if it was in the format '(0..X)', then that tells us how many
        % values there should be - we can make use of that as a check
        
        % we can either check this for every field, or check it for every
        % time we find a line with an empty field - arguably this is
        % simpler because the frequencies are of a similar order of
        % magnitude
        arrayParameters = sscanf(output{counter, 2}, '(%d..%d)')';
        
        if numel(arrayParameters == 2)
            % then we got a match
            % this covers zero and 1-indexing
            output{counter, 3} = sum(arrayParameters) + 1;
            
            % and we can overwrite it
            output{counter, 2} = parameters{m};
            
        else
            % otherwise append it - the last one was apparently a previous
            % array
            output{counter, 2} = [output{counter, 2}, ' ', parameters{m}];
        end
        
    else
        % then the field name is good - keep it
        counter = counter + 1;
        output{counter, 1} = parameters{m};
        output{counter, 2} = values{m};
        output{counter, 3} = 1;
    end
end

% cleanup the unused fields at the end
output(counter + 1:end, :) = [];


% now we can go back over the array and cleanup properly
for m = 1:size(output, 1)
    
    % strip off the leading $ that indicate "custom" fields
    output{m, 1} = output{m, 1}(2:end);
    
    % if the first character of a field is '<', it is a string
    if strncmp(output{m, 2}, '<', 1)
        
        % format we need to process is:
        % <string1> <string2> <string3>...
        % the strtok example works on HTML but gives the contents of the
        % angle brackets as well what's in between (which we don't want)
        
        % there is probably a very clever way of doing this using regexp
        % but I couldn't get it, or sscanf, or textscan to work they way I
        % wanted it to
        
        % this also strips off the angle brackets in the process
        startPoint = find(ismember(output{m, 2}, '<')) + 1;
        endPoint = find(ismember(output{m, 2}, '>')) - 1;
        startValues = numel(startPoint);
        
        if startValues ~= numel(endPoint)
            warning('Undefined behaviour for additional angle brackets.')
        end
        
        tempCell = cell(startValues, 1);
        
        for n = 1:startValues
            tempCell{n} = output{m, 2}(startPoint(n):endPoint(n));
        end
        
        % turn it into a character vector if there's just one of them
        if numel(tempCell) == 1
            tempCell = tempCell{1};
            
            % if that is an empty (1x0) char, then convert it into a
            % standard format
            if isempty(tempCell)
                tempCell = '';
            end
        end
        
        output{m, 2} = tempCell;
        
    else
        
        % convert numbers where possible
        % some of the new parameters (LOCKED, LOCSHIFT, NUSFPNZ,
        % SCALEDBYNS, SCALEDBYRG) have strings, but no angle brackets
        % ...otherwise we could just have used the same logic for
        % everything - they all say "no" at the moment - perhaps this is a
        % yes/no exception?
        if output{m, 3} == 1
            possibleNumber = str2double(output{m, 2});
            
            if ~isnan(possibleNumber)
                output{m ,2} = possibleNumber;
            end
            
        else
            % convert the arrays - these are better-behaved
            output{m, 2} = str2num(output{m, 2});
            
            % check it worked
            if numel(output{m, 2}) ~= output{m, 3}
                warning('Expected size of array for %s (%d) was not the expected size (%d).', output{m, 1}, numel(output{m, 2}), output{m, 3})
            end
        end
    end
end

% and finally remove the last column - we don't need it anymore
output(:, 3) = [];

% and convert that into a structure
output = cell2struct(output(:, 2), genvarname(output(:, 1)));

allData.Notes = output;

% ------------------------------------------------------------------------------
function output = jcampread2(filename,warnFlag)
%JCAMPREAD2 reads JCAMP-DX (DX) files
%
%   This is a slightly modified version of JCAMPREAD that is more relaxed
%   about data structures. Specifically it is more tolerant of the header
%   group, and will also import files even though they don't contain
%   XYDATA. This is the case with Bruker TopSpin files.
%
%   JCAMPDATA = JCAMPREAD2(FILE) reads in JCAMP-DX format data from FILE
%   and creates a structure JCAMPDATA, containing these fields:
%           Title
%           DataType
%           DataClass (version 5.00 and above, optional)
%           Origin
%           Owner
%           Blocks (optional)
%           Notes
%
%   The Blocks field of the structure is an array of structures
%   corresponding to each set of data in the file. These structures have
%   the following fields:
%           XData
%           YData
%           ZData (if multiple blocks)
%           XUnits
%           YUnits
%           ZUnits (if multiple blocks)
%           Notes
%
%   The function supports Versions 4.24 and 5 of the JCAMP-DX
%   format for infrared and mass spectrometry.
%
%   For more details on the JCAMP-DX format, see
%   http://www.jcamp-dx.org/index.html.
%
%   Sample data is available from http://www.jcamp-dx.org/testdata.html.
%
%   Example:
%
%       % Read in a sample JCAMP-DX file and plot the frequency spectrum.
%       jcampStruct = jcampread('brukaffn.dx')
%       data = jcampStruct.Blocks(1);
%       plot(data.XData,data.YData);
%       title(jcampStruct.Title);
%       xlabel(data.XUnits);
%       ylabel(data.YUnits);
%
%   See also MSLOWESS, MSSGOLAY, MSVIEWER, MZXMLREAD.

%  Note: http://www.jcamp-dx.org
%  Copyright 1997-2007  International Union of Pure and Applied Chemistry

% Copyright 2003-2017 The MathWorks, Inc.

% Check first if the input is a file identifier (used for the compound
% structure) or a filename string.
openedFile = 0;
if nargin < 2
    warnFlag = false;
end
try
    pos = ftell(filename);
catch 
    pos = -1;
end
if (pos < 0)
    try
        fid = fopen(filename,'rt');
        openedFile = 1;
    catch 
        fid = -1;
    end
    if fid == -1
        error(message('bioinfo:jcampread:CannotOpenJCAMPFile',filename));
    end
else
    fid = filename;
end
output = [];

%=== Read in first 6 lines and verify the headers are present
headerField = {'##TITLE','##JCAMP-DX','##DATA TYPE','##DATA CLASS','##ORIGIN','##OWNER','##NUM DIM'};
headerField2 = {'##TITLE','##JCAMPDX','##DATATYPE','##DATACLASS','##ORIGIN','##OWNER','##NUMDIM'};
allHeaders = union(headerField,headerField2);
headers = cell(size(headerField));
missingHeaderWarning = [];
count = 1;
while count <= 6
    % Retrieve the Labeled Data Record
    pos = ftell(fid);
    line = getnext(fid);
    [field,rem] = strtok(line,'='); 
    field = strtrim(field);
    if strcmp(field,'##')||~strncmp(field,'##',2)
        continue;
    elseif ~ismember(field,allHeaders)
        if warnFlag
            % Create the warning message but wait to throw it, as the headers might appear later
            missingHeaderWarning.Message = message('bioinfo:jcampread:BadJCAMPHeader', headerField{count}, field); 
            missingHeaderWarning.FieldId = count;
            warnFlag = false;
        end
        fseek(fid,pos,-1);
        break;
    end
    
    [found, position] = ismember(field,headerField);
    if ~found
        [~, position] = ismember(field,headerField2);
    end
    % Remove whitespaces and comments from the value
    value = strtok(rem,'=');
    headers{position} = strtrim(strtok(value,'$'));

    % Skip the Data Class field if the version is earlier than v.5
    % use str2double to deal with [].
    if (position==3) && (str2double(headers{2})<5)
        count = 4;
    end
    count = count+1;
end

%=== Check whether the version of JCAMP file format is supported
if ~isempty(headers{2})
    version = sscanf(headers{2},'%f');
    if (version < 5)&&(version ~= 4.24)
        warning(message('bioinfo:jcampread:BadJCAMPVersion', sprintf( '%0.2f', version )));
    end
    
%=== Skip unsupported format JCAMP-CS 
elseif strcmp(field,'##JCAMP-CS')
    
    % Search for the beginning of the next block
    while ~feof(fid) && ~strcmp(field,'##END')
        line = getnext(fid);
        field = strtok(line,'=');
        field = strtrim(field);
    end
    return
end

%=== Parse the rest of Labeled Data Records
numBlocks = 0;
DataSpec = struct('XFactor',1,'YFactor',1,'XUnits','','YUnits','',...
    'NPoints',-1,'FirstX',-1,'LastX',-1,'FirstY',-1);
nCounter = 1; % Keeps track of the notes index
notes = cell(1,2);
inNote = false;
while ~feof(fid)
    
    %=== Make sure it is a vald line and extract field/value pairs
    pos = ftell(fid);
    line = getnext(fid);

    try
        % Some files have lines with multiple '='. Set field to the
        % text before the first '=' and value as everything after that.
        fieldAndValue = regexp(line,'##([^=]*)=(.*)','tokens');
        field = fieldAndValue{1}{1};
        value = fieldAndValue{1}{2};
        inNote = false;
    catch 
        % some files have multi-line notes without the ## delimiter
        if inNote
            notes{nCounter,1} = line;
            nCounter = nCounter+1;
        end
        continue;
    end
    try
        field = upper(strtrim(field));
    catch 
        continue;
    end
    try
        value = strtrim(value);
    catch 
        value = '';
    end
    
    
    switch field
        case 'XUNITS'
            DataSpec.XUnits = value;
        case 'YUNITS'
            DataSpec.YUnits = value;
        case 'ORIGIN'
            if isempty(headers{5})
                headers{5} = value;
                % Avoid throwing header warning if header is found before the data
                if ~isempty(missingHeaderWarning) && (missingHeaderWarning.FieldId == 5)
                    missingHeaderWarning = [];
                end
            end
        case 'OWNER'
            if isempty(headers{6})
                headers{6} = value;
                % Avoid throwing header warning if header is found before the data
                if ~isempty(missingHeaderWarning) && (missingHeaderWarning.FieldId == 6)
                    missingHeaderWarning = [];
                end
            end
        case 'XFACTOR'
            DataSpec.XFactor = sscanf(value,'%f');
        case 'YFACTOR'
            DataSpec.YFactor = sscanf(value,'%f');
        case 'FIRSTX'
            DataSpec.FirstX = sscanf(value,'%f');
        case 'LASTX'
            DataSpec.LastX = sscanf(value,'%f');
        case 'NPOINTS'
            DataSpec.NPoints = sscanf(value,'%d');
        case 'FIRSTY'
            DataSpec.FirstY = sscanf(value,'%f');
        case 'BLOCKS'
            numBlocks = sscanf(value,'%f');
        case 'TITLE'
            % Run subroutine for compound data
            if ~isempty(missingHeaderWarning)
                warning(missingHeaderWarning.Message); % Reached a new block, throw header warning
                missingHeaderWarning = [];
            end
            fseek(fid, pos, 'bof'); % Move back one line to start at TITLE
            
            % LINK block data may be unsupported, like a JCAMP-CS data block
            linkFlag = strcmp(headers{3},'LINK');
            dataBlocks = compoundRead(fid, numBlocks, linkFlag, warnFlag);
            break;
        case 'XYDATA'
            % Run subroutine for xydata
            if ~isempty(missingHeaderWarning)
                warning(missingHeaderWarning.Message); % Reached data, throw header warning
                missingHeaderWarning = [];
            end
            dataBlocks(1) = xyRead(fid,DataSpec);
            break;
        case'NTUPLES'
            % Run subroutine for ntuples
            if ~isempty(missingHeaderWarning)
                warning(missingHeaderWarning.Message); % Reached data, throw header warning
                missingHeaderWarning = [];
            end
            numBlocks = 1; %#ok Multiple dataBlocks, so expect an additional END field
            dataBlocks = ntupleRead(fid);
            break;
        case 'PEAK TABLE'
            % Run subroutine for peak table
            if ~isempty(missingHeaderWarning)
                warning(missingHeaderWarning.Message); % Reached data, throw header warning
                missingHeaderWarning = [];
            end
            dataBlocks(1) = peakRead(fid,DataSpec);
            break;
        case 'XYPOINTS'
            % Run subroutine for xypoints
            if ~isempty(missingHeaderWarning)
                warning(missingHeaderWarning.Message); % Reached data, throw header warning
                missingHeaderWarning = [];
            end
            dataBlocks(1) = xyPointsRead(fid,DataSpec);
            break;
%         case 'END'
%             % If you find an END without a multiple dataBlock structure, there
%             % is an error in the file format.
%             if numBlocks < 1
%                 % About to error out, throw header warning first
%                 if ~isempty(missingHeaderWarning)
%                     warning(missingHeaderWarning.Message);
%                 end
%                 error(message('bioinfo:jcampread:BadJCAMPDataType', headers{ 3 }))
%             end
        otherwise
            % Store field/value as a note
            notes{nCounter,1}=field;
            notes{nCounter,2}=value;
            nCounter = nCounter+1;
            % some example files have multi-line notes without the ##
            % delimiter
            inNote = true;
    end

end

if ~isempty(missingHeaderWarning)
    warning(missingHeaderWarning.Message); % No data was reached, throw header warning first
end

if(openedFile)
    fclose(fid);
end

% if ~exist('dataBlocks','var')
%     error(message('bioinfo:jcampread:BadJCAMPDataType', headers{ 3 }))
% end

output.Title = headers{1};
output.DataType = headers{3};
if ~isempty(headers{4})
    output.DataClass = headers{4};
end
output.Origin = headers{5};
output.Owner = headers{6};

if exist('dataBlocks','var')
  output.Blocks = dataBlocks;
end

output.Notes = notes;

% -----------------------------------------------------------------
% PEAKREAD
% Handles uncompressed data in the (XY..XY) format
% -----------------------------------------------------------------
function dataBlock = peakRead(fid,DataSpec)
if (DataSpec.NPoints>0)
    x = zeros(1,DataSpec.NPoints);
    y = zeros(1,DataSpec.NPoints);
else
    x = [];
    y = [];
end

j = 1;
pos = ftell(fid);
data = getnext(fid);
% some files have delimiter ',' so we need to be robust to that.
formatString = '%f, %f';
delimiter = ';';
try % try the default format of commas and semi-colons.
    textscan(data,formatString,'delimiter',delimiter);
catch allExceptions %#ok<NASGU>
    try % try using comma as delimiter for everything
        formatString = '%f %f';
        delimiter = ',';
        textscan(data,formatString,'delimiter',delimiter);
    catch  % go back to default
        warning(message('bioinfo:jcampread:BadDelimiter'));
        formatString = '%f, %f';
        delimiter = ';';
    end
end
while ~strncmp(data,'##',2)
    xy = textscan(data,formatString,'delimiter',delimiter);
    xi = xy{1};
    yi = xy{2};
    
    k = j+length(xi)-1;
    x(j:k) = xi;
    y(j:k) = yi;
    j = k+1;
    pos = ftell(fid);
    data = getnext(fid);
end
% Don't skip the first line of the next block, if it exists
if ~strcmp(strtok(data,'='),'##END')
    fseek(fid,pos,'bof');
end

dataBlock.XData = x*DataSpec.XFactor;
dataBlock.YData = y*DataSpec.YFactor;
dataBlock.XUnits = DataSpec.XUnits;
dataBlock.YUnits = DataSpec.YUnits;

% -----------------------------------------------------------------
% XYREAD
% Handles data in the (X++(Y..Y)) form, including SQZ, DIF, and DUP
% compressed data.
% -----------------------------------------------------------------
function dataBlock = xyRead(fid,DataSpec)
tol = 1; % Set the tolerance for the checkpoints
xTolWarnFlag = true;
yTolWarnFlag = true;

SQZ = 'ihgfedcba@ABCDEFGHI';
DIF = 'rqponmlkj%JKLMNOPQR';
DUP = 'STUVWXYZs';

if (DataSpec.NPoints>0)
    x = linspace(DataSpec.FirstX/DataSpec.XFactor,...
        DataSpec.LastX/DataSpec.XFactor,DataSpec.NPoints);
    y = zeros(1,DataSpec.NPoints);
else
    error(message('bioinfo:jcampread:BadJCAMPDataPoints'));
end

j=1;
mode = 1;

pos = ftell(fid);
data = getnext(fid);
while ~strncmp(data,'##',2)
  values = regexp(data,'[+-]*[a-s_A-Z_%@]?\d*\.?\d*([Ee][+-]\d+)?','match');
    xval = sscanf(values{1},'%f');
    if abs(xval-round(x(j)))>tol
        if xTolWarnFlag
            warning(message('bioinfo:jcampread:XOutofTol'))
            xTolWarnFlag = false;
        end
    end
    for i = 2:length(values)
        num = values{i};
        if contains(SQZ,num(1)) % ~isempty(strfind(SQZ,num(1)))
            dig = find(SQZ==num(1))-10;
            num = sscanf([sprintf('%d',dig) num(2:end)],'%f');
            if (i==2)&&(j>1)&&(mode==2)
                if abs(num-y(j))>tol
                    if yTolWarnFlag
                        warning(message('bioinfo:jcampread:YOutofTol'))
                        yTolWarnFlag = false;
                    end
                    y(j) = num;
                end
            else
                y(j) = num;
            end
            j = j+1;
        elseif contains(DIF,num(1)) % ~isempty(strfind(DIF,num(1)))
            difDig = find(DIF==num(1))-10;
            difValue = sscanf([sprintf('%d',difDig) num(2:end)],'%f');
            y(j) = y(j-1)+difValue;
            j = j+1;
            mode = 2;
        elseif contains(DUP, num(1)) %~isempty(strfind(DUP,num(1)))
            dupDig = find(DUP==num(1));
            if length(num)==1
                dupValue = dupDig;
            else
                dupValue = sscanf([sprintf('%d',dupDig) num(2:end)],'%f');
            end
            for k = 1:dupValue-1
                if mode==1
                    y(j) = y(j-1);
                else
                    y(j) = y(j-1)+difValue;
                end
                j = j+1;
            end
        else
            if (i==2)&&(j>1)&&(mode==2)
                if abs(sscanf(num,'%f')-y(j))>tol
                    warning(message('bioinfo:jcampread:YOutofTol'))
                end
            else
                y(j) = sscanf(num,'%f');
            end
            j = j+1;
        end
    end
    if (mode==2)
        j = j-1; % rewind one value for the checkpoint
    end
    pos = ftell(fid);
    data=getnext(fid);
end
% Don't skip the first line of the next block, if it exists
if ~strcmp(strtok(data,'='),'##END')
    fseek(fid,pos,'bof');
end

dataBlock.XData = x*DataSpec.XFactor;
dataBlock.YData = y*DataSpec.YFactor;
dataBlock.XUnits = DataSpec.XUnits;
dataBlock.YUnits = DataSpec.YUnits;

%------------------------------------------------------------------
% XYPOINTSREAD
% Reads XY pairs of spectral data not spaced at constant abscissa
% increments. This is intended for infrared filter spectrometers.
%------------------------------------------------------------------
function dataBlock = xyPointsRead(fid,DataSpec)
if (DataSpec.NPoints>0)
    x = zeros(1,DataSpec.NPoints);
    y = zeros(1,DataSpec.NPoints);
else
    x = [];
    y = [];
end

j = 1;
pos = ftell(fid);
data = getnext(fid);

% X and Y are separated by commas; XY pairs are separated by semicolons
% or blanks.
formatString = '%f, %f'; 
delimiter = ';';
try % try the default format of commas and semi-colons.
    textscan(data,formatString,'delimiter',delimiter);
catch 
    try % try using blank as delimiter
        delimiter = ',';
        textscan(data,formatString,'delimiter',delimiter);
    catch % go back to default
        warning(message('bioinfo:jcampread:BadDelimiter'));
        delimiter = ';';
    end
end
while ~strncmp(data,'##',2)
    xy = textscan(data,formatString,'delimiter',delimiter);
    xi = xy{1};
    yi = xy{2};
    
    k = j+length(xi)-1;
    x(j:k) = xi;
    y(j:k) = yi;
    j = k+1;
    pos = ftell(fid);
    data = getnext(fid);
end
% Don't skip the first line of the next block, if it exists
if ~strcmp(strtok(data,'='),'##END')
    fseek(fid,pos,'bof');
end

dataBlock.XData = x*DataSpec.XFactor;
dataBlock.YData = y*DataSpec.YFactor;
dataBlock.XUnits = DataSpec.XUnits;
dataBlock.YUnits = DataSpec.YUnits;

% -----------------------------------------------------------------
% COMPOUNDREAD
% This subfunction recursively call the main function to read in each
% subblock of data.
% -----------------------------------------------------------------
function dataBlocks = compoundRead(fid,numBlocks,linkFlag,warnFlag)
dataBlocks = [];
if nargin < 4
    warnFlag = true;
    if nargin < 3
        linkFlag = false;
    end
end
for i = 1:numBlocks
    try
        DataStr = jcampread2(fid,warnFlag);
        if linkFlag && isempty(DataStr)
            % LINK block may have been unsupported JCAMP-CS format
            continue
        end
        try
            dataBlocks(i) = DataStr.Blocks(1); %#ok<*AGROW>
        catch allExceptions %#ok<NASGU>
            dataBlocks(i).XData = DataStr.Blocks(1).XData;
            dataBlocks(i).YData = DataStr.Blocks(1).YData;
            dataBlocks(i).XUnits = DataStr.Blocks(1).XUnits;
            dataBlocks(i).YUnits = DataStr.Blocks(1).YUnits;
        end
        dataBlocks(i).Title = DataStr.Title;
        dataBlocks(i).DataType = DataStr.DataType;
        dataBlocks(i).Owner = DataStr.Owner;
        dataBlocks(i).Origin = DataStr.Origin;
        dataBlocks(i).Notes = DataStr.Notes;
    catch theException
        msgstr = theException.message;
        if warnFlag
            warning(message('bioinfo:jcampread:CompoundStructureError', msgstr));
            warnFlag = false;
        end
        while (~feof(fid))&&(~strncmp(fgetl(fid),'##END',5))
        end
    end
    if ~exist('dataBlocks','var')
        dataBlocks = struct;
    end
end

% -----------------------------------------------------------------
% NTUPLEREAD
% Handles Spectral Series Data
% -----------------------------------------------------------------
function dataBlocks = ntupleRead(fid)
blockNum = 1;
DataSpec = struct('XFactor',1,'YFactor',1,'XUnits','','YUnits','',...
    'NPoints',-1,'FirstX',-1,'LastX',-1,'FirstY',-1);
data = getnext(fid);
while ~strncmp(data,'##END',5)
    try
        t=regexp(data,'##([\s\w-/]*)\s*=\s*(.*)','tokens');
        field = t{1}{1};
        value = t{1}{2};
    catch 
        data = getnext(fid);
        continue;
    end

    if strncmp(field,'$',1)
        continue; % Ignore custom fields designated by ##$...
    elseif strcmp(field,'VAR_NAME')
        name = textscan(value,'%s','delimiter',',');
        name = name{1};
    elseif strcmp(field,'SYMBOL')
        cellsymbol = textscan(value,'%s','delimiter',',');
        cellsymbol = cellsymbol{1};
    elseif strcmp(field,'VAR_DIM')
        vardim = textscan(value,'%d','delimiter',',');
        vardim = vardim{1};
    elseif strcmp(field,'UNITS')
        units = textscan(value,'%s','delimiter',',');
        units = units{1};
    elseif strcmp(field,'FIRST')
        first = textscan(value,'%f','delimiter',',');
        first = first{1};
    elseif strcmp(field,'LAST')
        last = textscan(value,'%f','delimiter',',');
        last = last{1};
    elseif strcmp(field,'FACTOR')
        fctrs = textscan(value,'%f','delimiter',',');
        fctrs = fctrs{1};
    elseif strcmp(field,'PAGE')
        z = textscan(value,'%s %f','delimiter','=');
        zvar = strtrim(z{1}{1}); % Extract one character
        zval = z{2};
    elseif strcmp(field,'NPOINTS')
        DataSpec.NPoints = sscanf(value,'%d');
    elseif strcmp(field,'DATA TABLE')
        
        % Determine the type of data
        form = textscan(value,'%s','delimiter',',');
        form = form{1};
        mode = 0;
        vars = regexp(form{1},'\((\w)(\w)..\1\2\)','tokens');
        if isempty(vars)
            mode = 1;
            vars = regexp(form{1},'\((\w.*)\+\+\((\w)..\2\)','tokens');
            if isempty(vars)
                error(message('bioinfo:jcampread:BadNtupleFormat', form{ 1 }));
            end
        end
        if ~isempty(cellsymbol)
            Xidx = find(strcmp(cellsymbol,vars{1}{1}));
            Yidx = find(strcmp(cellsymbol,vars{1}{2}));
            Zidx = find(strcmp(cellsymbol,zvar));
        else
            error(message('bioinfo:jcampread:BadJCAMPntuple'));
        end
        if exist('fctrs','var')&&~isempty(fctrs)
            DataSpec.XFactor = fctrs(Xidx);
            DataSpec.YFactor = fctrs(Yidx);
        end
        if exist('units','var')&&~isempty(units)
            DataSpec.XUnits = units{Xidx};
            DataSpec.YUnits = units{Yidx};
        end
        if exist('first','var')&&~isempty(units)
            DataSpec.FirstX = first(Xidx);
            DataSpec.FirstY = first(Yidx);
        end
        if exist('last','var')&&~isempty(last)
            DataSpec.LastX = last(Xidx);
        end
        if DataSpec.NPoints < 0
            try
                DataSpec.NPoints = vardim(Xidx);
            catch allExceptions
                error(message('bioinfo:jcampread:BadJCAMPDataPoints'));
            end
        end
        
        % Run subroutine based on table format
        if mode
            dataBlock = xyRead(fid,DataSpec);
        else
            dataBlock = peakRead(fid,DataSpec);
        end
        
        % Restructure result and add extra values
        if exist('name','var')&&~isempty(name)
            dataBlock.XName = name{Xidx};
            dataBlock.YName = name{Yidx};
            dataBlock.ZName = name{Zidx};
        end
        
        try
            dataBlock.ZUnits = units{Zidx};
        catch allExceptions %#ok<NASGU>
            dataBlock.ZUnits = '';
        end
        dataBlock.ZData = zval;
        dataBlocks(blockNum) = dataBlock;
        blockNum = blockNum+1;
        
        % Move back up to the end of the table so next line is Page or END
%         fseek(fid,endTable,'bof');
    end
    data = getnext(fid);
end

% -----------------------------------------------------------------
% GETNEXT - gets the next valid line from the file,
% ignores comments and custom fields.
% -----------------------------------------------------------------
function data = getnext(fid)
data = fgetl(fid);
%while strncmp(data,'$$',2)||strncmp(data,'##$',3)||isempty(data)
while strncmp(data,'$$',2)||isempty(data)
    data = fgetl(fid);
end
data = strtrim(data);
if data<0
    error(message('bioinfo:jcampread:NoJCAMPEnd'))
% else
%     if contains(data,'$') % ~isempty(strfind(data,'$'))
%         data = textscan(data,'%s',1,'delimiter','$');
%         data = strtrim(data{1});
%         data = data{1};
%     end
end

function ret=contains(A,B)

  ret = ~isempty(strfind(A,B));
  
function s=message(varargin)
  s = sprintf('%s ', varargin{:});
        
