function data = read_mccode(filename)
% read_mccode: load a McCode (McStas/McXtrace) simulation result
%   data = read_mccode(filename)
%
% This function imports a McCode simulation result
% It returns the McStas simulation structures. 
% It may import both McCode and Matlab format data sets.
%
% input:
%  filename: one or more simulation name(s) or directory
%          or a single detector file name
%          if filename does not exist, a file selector is called.
% 
% output:
%  data: a cell of detector structures
%
% examples:
%   read_mccode mcstas.sim
%
% (c) E.Farhi, ILL. License: EUPL.
% See also: read_anytext, read_inx


% Written by: E. Farhi
% Date: April 16st 2010
% Release: McStas 1.6
% Version: $Revision$
% Origin: ILL
%
%   This file is part of the McCode neutron/X ray-trace simulation package
%   Copyright (C) 1997-2010, All rights reserved
%   Risoe National Laborartory, Roskilde, Denmark
%   Institut Laue Langevin, Grenoble, France
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; version 2 of the License.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%

% check for input argument: filename ?
  data = {};
  
  if nargin == 0
    return
  end

  % import data set
  data = mcplot_load(filename);
  
  % check data structures
  data = mcplot_check_data(data);

end % read_mccode (main)

% ==============================================================================
% inline functions
% ==============================================================================

% function match= mcplot_filestrfind(filename, pattern, {buffer})
% function match= mcplot_filefgetl(filename, positions)
% function data = mcplot_load(filename)
% function data = mcplot_load_mccode(filename)
% function data = mcplot_load_matlab(filename)
% function data = mcplot_load_structure(s)
% function data = mcplot_load_sim(filename)
% function data = mcplot_load_scan(filename)
% function data = mcplot_split_multiarray(structure)
% function mcplot_display(structs) 

% ==============================================================================

function match = mcplot_filestrfind(filename, string, buffer)
  % mcplot_filestrfind: find matches of 'string' in the specified filename
  % filename: name of file to search into
  %   string: pattern to search for
  %   buffer: size of buffer used to parse file. buffer=0 will read the entire file
  %           a finite buffer size will only read that amount within the file,
  %           corresponding e.g. with the header.
  %   RETURN: positions of matches within the file
  
  match = [];
  
  if nargin < 2, return; end
  if nargin < 3, 
    buffer = Inf; % will read whole file iteratively
  end
  if buffer == Inf, buffer   = 0;      readToEOF= 1; else readToEOF = 0; end
  if ~buffer,       buffer   = 10000;  end  % default buffer size
  
  invalid         = find(~isstrprop(string,'print'));
  string(invalid) = ' ';  % replace non printable characters with spaces
  
  if iscellstr(filename), filename = filename{1}; end
  
  fid = fopen(filename);
  if fid == -1, return; end   % also returns when filename is empty
  
  offset = 0; % offset at which the buffer must be loaded from file
  
  % read file by blocks, once or iteratively
  while (~feof(fid) && isempty(ferror(fid))) 
    block = fread(fid, buffer, 'uint8=>char');
    block = block(:)';
    if ~isempty(block)
      invalid        = find(~isstrprop(block, 'print'));
      block(invalid) = ' '; % replace non printable characters with spaces
      
      % find search string in the block
      this_match = strfind(block, string);      % search string
      if ~isempty(this_match)
        match = [ match (this_match+offset) ];  % append new match locations (absolute)
      end
      
      % prepare new offset for iterative search
      offset = ftell(fid);
      if offset == -1, readToEOF = 0; end
    else 
      readToEOF = 0; 
    end
    if readToEOF == 0, break; end
  end % while (readToEOF)
  
  fclose(fid);  % close file
  
end % mcplot_filestrfind

% ==============================================================================

function match = mcplot_filefgetl(filename, positions)
  % mcplot_filefgetl: extract iteratively lines at given positions in file
  %  filename: name of file to search into
  % positions: positions of matches within the file
  % RETURN: cellstr of lines corresponding to positions in file
  
  match = {};
  
  fid = fopen(filename);
  if fid == -1, return; end % also returns when filename is empty
  
  match = cell(1,length(positions));
  
  for index=1:length(positions)
    fseek(fid, positions(index), -1);  % jump at match position
    line = fgetl(fid);
    if ~ischar(line), break; end
    match{index} = strtrim(line);
  end
  
  fclose(fid);  % close file
end % mcplot_filefgetl

% ==============================================================================

function data=mcplot_load(filename)
  % mcplot_load: loads data from directory or file name
  %   Handles McCode and Matlab format, as overview (.sim), directory and single data file
  % filename: name of file to load, or a structure
  %   RETURN: a single structure or a cell array of data structures
  
  data = {};
  
  if ishandle(filename)
    data = get(filename, 'UserData');
  elseif isstruct(filename)
    [data,parameters] = mcplot_load_structure(filename);
  elseif iscellstr(filename)
    for index=1:length(filename)
      data = { data{:} mcplot_load(filename{index}) };
    end
    return
  elseif isdir(filename)  % search for mcstas.sim, mcstas.dat and content.sim in directory
    if ~isempty(     dir(fullfile(filename, 'mcstas.sim'))) % original mcstas files
      data = mcplot_load(fullfile(filename, 'mcstas.sim'));
    elseif ~isempty( dir(fullfile(filename, 'mcstas.dat')))
      data = mcplot_load(fullfile(filename, 'mcstas.dat'));
    elseif ~isempty( dir(fullfile(filename, 'mccode.sim'))) % common mcstas/mcxtrace files
      data = mcplot_load(fullfile(filename, 'mccode.sim'));
    elseif ~isempty( dir(fullfile(filename, 'mccode.dat')))
      data = mcplot_load(fullfile(filename, 'mccode.dat'));
    elseif ~isempty( dir(fullfile(filename, 'mcstas.m')))   % matlab related files
      data = mcplot_load(fullfile(filename, 'mcstas.m'));
    elseif ~isempty( dir(fullfile(filename, 'mccode.m')))
      data = mcplot_load(fullfile(filename, 'mccode.m'));
    end
    % else get all files in directory
  end
  if isempty(data)
  
    this_dir = dir(filename);           % get file reference (and test for existence)
    if isempty(this_dir), return; end   % invalid. exit. Also returns when filename is empty
    
    filepath =fileparts(filename);      % determine directory
    if isempty(filepath), filepath = pwd; end
    
    % find directory entries that are not directories, and non NULL
    index    = find(real([this_dir.isdir]) == 0 & real([this_dir.bytes]) > 0);
    this_dir = char(this_dir.name);
    this_dir = (this_dir(index,:));
    rdir     = cellstr(this_dir);
    rdir     = strcat([ filepath filesep ], char(rdir));
    filename = cellstr(rdir);           % directory listing as cell (names only)

    % filename is now a cellstr: IMPORT HERE
    for index=1:length(filename)
      % test if this is a sim (contains 'begin instrument') or multiarray_1d file
      % in this case, the returned value is a full data set
      this_data = mcplot_load_sim(filename{index});       % load McCode SIM (overview) files
      if isempty(this_data)
        this_data = mcplot_load_scan(filename{index});    % load McCode scan files (multiarray_1d)
      end
      if isempty(this_data)
        this_data = mcplot_load_matlab(filename{index});  % load Matlab format files
      end
      if isempty(this_data)
        this_data = mcplot_load_mccode(filename{index});    % load single McCode monitor files
      end
      if isstruct(this_data) && ~isfield(this_data, 'filename') setfield(this_data,'filename',filename{index}); end
      data = { data{:} this_data };
    end
    
  end
  % in case result is embedded inside a cell, we extract it
  if iscell(data) && iscell(data{1}) && length(data) == 1
    data = data{1};
  end
  if iscell(data) && length(data) == 1
    data = data{1};
  end

end % mcplot_load

% ==============================================================================

function structure=mcplot_load_mccode(filename)
  % mcplot_load_mccode: load a single data set (not sim) or a structure
  % filename: name of file to load
  %   RETURN: a cell array of data structures
  
  structure = [];
  
  isMcCode = [ mcplot_filestrfind(filename, 'array_1d', 10000) ...
               mcplot_filestrfind(filename, 'array_2d', 10000) ...
               mcplot_filestrfind(filename, 'multiarray_1d', 10000) ];
  
  if isempty(isMcCode), return; end % also returns when filename is empty
  
  disp([ 'Loading ' filename ' (McCode format)' ]);
  % with scilab, a call to fscanfMat will extract data and header
  fid = fopen(filename);
  if fid == -1, return; end % also returns when filename is empty
  
  % read header
  header = textscan(fid,'#%s','endOfLine','\n','delimiter','\n');

  if iscellstr(header{1}) && length(header)==1
    header = header{1};
  end
  
  paramstr = '';
  % build structure from header fields
  for index=1:length(header)
    % use strtok to split line around':'
    [field, value] = strtok(header{index}, ':');
    value(find(value == ':')) = ''; value=strtrim(value); field=strtrim(field);
    value=strrep(value,'''','');  % remove quotes
    % keep 'value' as a string for further use before converting to num
    value_str=value;
    if all(ismember(value, sprintf(' 0123456789+-.eEdDi\n\r\t\f')))
      num   = str2num(value);
    else num=[]; end
    if ~isempty(num), value = num; end
    if strncmp(field,'Instrument', length('Instrument'))
      field = 'Instrument';
    else
      field = genvarname(field);                        % validate variable name
    end
    if isempty(strfind(field, 'Param'))
      if ~isempty(value), structure = setfield(structure, field, value); end   % set new field
    else  % special case for parameters. Build 'Param' sub structure
      if ~isfield(structure,'Param'), structure = setfield(structure, 'Param',[]); end  % create if needed
      param        = getfield(structure, 'Param');
      % use strtok to split line around'='
      [var,value]  = strtok(value_str,'=');
      value(find(value == '=')) = ''; value=strtrim(value); var=strtrim(var);
      paramstr = [ paramstr var '=' value '; ' ];
      num          = str2num(value);
      if ~isempty(num), value = num; end
      param        = setfield(param, var, value);
      structure    = setfield(structure, 'Param', param);     % store updated sub-structure
    end
  end
  structure = setfield(structure, 'Parameters', paramstr);
  structure.filename = filename;
  clear header
  
  % additional compatibility checks (from McCode/Matlab format)
  if isfield(structure, 'parameters') && ~isfield(structure,'Param')
    param = getfield(structure,'parameters');
    structure = rmfield(structure,'parameters');
    structure = setfield(structure, 'Param', param);
  end
  if ~isfield(structure,'Param')
    structure = setfield(structure, 'Param', 'Unknown instrument parameters');
  end
  
  % extract type and theoretical dimension of data block
  if isfield(structure,'type')
    % use strtok to split line around'()'
    [t s]=strtok(structure.type, '()');
    structure.type = t;
    % remove parenthesis
    s=strrep(s,'(',''); s=strrep(s,')','');
    structure.size = str2num(s);
  end
  
  % load data block
  if exist('textscan') && ~exist ('OCTAVE_VERSION', 'builtin')
    frewind(fid);
    data     = textscan(fid,'%f','CommentStyle','#');
  end % else was obtained in one call to mcplot_textscan
  
  fclose(fid);
  if iscell(data) && length(data)==1
    data = data{1};
  end
  
  structure.data = data;
  clear data
  
  % reshape data block
  if ~isempty(strmatch('array_1d',structure.type))
    if prod(size(structure.data))/prod(structure.size) == 4
      % textscan provides a single vector: we must reshape as 4 columns
      structure.data = transpose(reshape(structure.data,[ 4 structure.size ]));
      % extract signal, errors and events
      if size(structure.data,2) >= 4, structure.events=structure.data(:,4); end
      if size(structure.data,2) >= 3, structure.errors=structure.data(:,3); end
      structure.x     =structure.data(:,1);
      structure.data  =structure.data(:,2);
    end
  elseif ~isempty(strmatch('array_2d',structure.type)) % 2d
    % update size from the actual read data (can be different for e.g. event list)
    if rem(numel(structure.data), structure.size(1)) == 0 && numel(structure.data) ~= prod(structure.size)*3
      structure.size = [ structure.size(1) numel(structure.data)/structure.size(1) ];
    end
    len      = prod(structure.size);
    structure.errors=[];
    structure.events=[];
    if exist('textscan') && ~exist ('OCTAVE_VERSION', 'builtin')  % Matlab read file and produces a single vector to be re-organized
      this_data        = structure.data;
      structure.data   = transpose(reshape(this_data(1:len),structure.size));
      if prod(size(this_data)) >= 2*len
        structure.errors = transpose(reshape(this_data((len+1):(2*len)),structure.size));
      end
      if prod(size(this_data)) >= 3*len
        structure.events = transpose(reshape(this_data((2*len+1):(3*len)),structure.size));
      end
    else
      % we used our own data reader (slower), which already shaped data block
      l = structure.size(1);
      this_data        = structure.data;
      structure.data = this_data(1:l,1:structure.size(2));
      
      if prod(size(structure.data)) >= 2*len
        structure.errors = this_data((l+1):(2*l),1:structure.size(2));
      end
      if prod(size(structure.data)) >= 3*len
        structure.events = this_data((2*l+1):(3*l),1:structure.size(2));
      end
    end
  elseif ~isempty(strmatch('multiarray_1d',structure.type)) % multiarray_1d
    % reshape data block according to the multiarray_1d(dim)
    n = prod(size(structure.data))/prod(structure.size);
    structure.size= [ structure.size n ];
    structure.data = transpose(reshape(structure.data,fliplr(structure.size)));
  end

end % mcplot_load_mccode

% ==============================================================================

function data = mcplot_check_data(structure)
  % mcplot_check_data: check all data sets for consistency
  data={};
  if iscell(structure)
    for index=1:length(structure)
      data = { data{:} mcplot_check_data(structure{index}) };
    end
    return
  end
  
  % check a single structure format fields
  if isfield(structure, 'Source') && isempty(dir(structure.Source))
    % not a file name
    if ~isfield(structure, 'Instrument')
      structure.Instrument = structure.Source;
    end
  end

  if ~isfield(structure, 'component') structure.component = 'unknown'; end
  if ~isfield(structure, 'filename')  structure.filename  = pwd; end
  if ~isfield(structure, 'errors')    structure.errors    = [];  end
  if ~isfield(structure, 'events')    structure.events    = [];  end
  if ~isfield(structure, 'title'),    structure.title     = '';  end
  structure.Source = structure.filename;

  % extract type and theoretical dimension of data block
  if ~isfield(structure, 'size')
    structure.size = [];
  end
  if isfield(structure,'type')
    % use strtok to split line around'()'
    [t s]=strtok(structure.type, '()');
    structure.type = t;
    % remove parenthesis
    s=strrep(s,'(',''); s=strrep(s,')','');
    structure.size = str2num(s);
  end
  
  if isempty(structure.size)
    if isfield(structure,'data')
      structure.size = size(structure.data);
    else
      structure.size = [0 0];
    end
  end

  % reshape data block from 'type'
  
  data = structure;
  
  % dimensions
end % mcplot_check_data

% ==============================================================================

function data=mcplot_load_matlab(filename)
  % mcplot_load_matlab: test and open a McCode/Matlab data file and load its content
  % filename: name of file to load (matlab format)
  %   RETURN: a data structure, or empty if not a McCode/Matlab file
  
  data = {};
  
  % search for 'Matlab' in header. Return if not found
  isMatlabScript = [ mcplot_filestrfind(filename, 'Matlab', 10000) ...
                     mcplot_filestrfind(filename, 'matlab', 10000) ...
                     mcplot_filestrfind(filename, '% Embedded function for building', 10000) ];
  
  if isempty(isMatlabScript), return; end % also returns when filename is empty
  
  % in principle, McCode/Matlab data format contains its own 'mcload_inline' routine
  % make sure we get a Matlab .m file
  cur_dir = pwd;
  [pathname, object, ext]=fileparts(filename);
  if isempty(pathname), pathname=pwd; end
  object = genvarname(object);
  % create a local copy with right extension (removed afterwards)
  if ~strcmp(filename, fullfile(pathname, [ object '.m' ]))
    copyfile(filename, fullfile(pathname, [ object '.m' ]))
  end
  cd(pathname);  % go in directory where data is
  
  % evaluate Matlab format data file
  data = eval(object, '[]');
  if ~strcmp(filename, fullfile(pathname, [ object '.m' ]))
    delete([ object '.m' ]);
  end
  cd(cur_dir);
  data.filename = filename;
  disp([ 'Loading ' filename ' (Matlab format)' ]);

  [data, parameters] = mcplot_load_structure(data);   % extract monitors recursively

  % insert extracted parameteres in each monitor structure
  if isstruct(data), 
    data.Param=parameters;
  elseif iscell(data)
    for index=1:length(data)
      this_data = data{index};
      if isstruct(this_data), 
        this_data.Param=parameters; 
        data{index} = this_data; 
      end
    end
  end
  
end % mcplot_load_matlab

% ==============================================================================

function [data,parameters]=mcplot_load_structure(s,parameters)
  % mcplot_load_structure: load data from a structure recursively
  %         s: structure with possibly full hierachy. We then search for those that contain some 'data'
  %    RETURN: a single structure or a cell array of data structures
  
  if nargin<2,
    parameters=[];
  end
  data = {};
  
  if ~isstruct(s), return; end
  
  if isfield(s,'data') && isnumeric(s.data)                  % found 'data': we keep that
    data = s;
    data.Param = parameters;
    disp([ 'Loading ' data.filename ]);
    data = { data };
    % return it as a structure
  elseif isfield(s,'name') && strcmp(s.name, 'parameters')
    parameters = s;
  else
    tag_names = fieldnames(s);
    for index=1:length(tag_names)       % scan recursively structure fields
      d = getfield(s,tag_names{index});
      if isstruct(d)
        [this_data,parameters] = mcplot_load_structure(d,parameters);
        if ~isempty(this_data)
          if iscell(this_data) && length(this_data) == 1 && isstruct(this_data{1})
            this_data = this_data(1);
          end
          data = { data{:} this_data{:} };
        end
      end
    end
  end
  if iscell(data) && length(data) == 1 && isstruct(data{1}) 
    data=data(1);
  end
end % mcplot_load_structure
  
% ==============================================================================

function data=mcplot_load_sim(filename)
  % mcplot_load_sim: test and open a .sim (overview) McCode data file and load its content
  % filename: name of file to load
  %   RETURN: a cell array of data structures, or empty if not a scan file
  
  data = {};
  
  % search for 'begin instrument' in header. Return if not found
  isSimFile    = mcplot_filestrfind(filename, 'begin instrument', 10000);
  if isempty(isSimFile)
    return  % also returns when filename is empty
  end
  disp([ 'Loading ' filename ' (McCode simulation overview)' ]);
  % search for 'filename:' tags, and extract links to other files
  filenameReferences = mcplot_filestrfind(filename, 'filename:');
  filenameLines      = mcplot_filefgetl  (filename, filenameReferences+length('filename:'));
  filepath = fileparts(filename);

  % loop on filenames
  for index = 1:length(filenameLines)
    % calls mcplot_load_mccode
    this_filename = fullfile(filepath, filenameLines{index});
    this_data     = mcplot_load_mccode(this_filename);
    if ~isempty(this_data) && ~isempty(strmatch('multiarray_1d',this_data.type))
      this_data.size = size(this_data.data);
      this_data = mcplot_split_multiarray(this_data);
    end
    data = { data{:} this_data }; % IMPORT sim file references
  end

end % mcplot_load_sim

% ==============================================================================

function data=mcplot_load_scan(filename)
  % mcplot_load_scan: test and open a multiarray (scan) McCode data file and load its content
  % filename: name of file to load
  %   RETURN: a cell array of data structures, or empty if not a sim/scan file
  
  data = {};
  
  % search for 'multiarray_1d' in header. Return if not found
  isMultiArray = mcplot_filestrfind(filename, 'multiarray_1d',    10000);
  if isempty(isMultiArray)  
    return  % also returns when filename is empty
  end
  
  % % IMPORT multiarray file
  data = mcplot_load_mccode(filename);                

  % then extract columns in a loop
  data = mcplot_split_multiarray(data);
  
end % mcplot_load_scan

% ==============================================================================

function data = mcplot_split_multiarray(structure)
  % mcplot_split_multiarray: load a multiarray data set and generate monitor column structures
  % structure: single data set structure (multiarray type)
  %    RETURN: a cell array of data structures, one for each monitor in the multiarray
  
  data = {};
  % first check if this is a multiarray
  if isempty(strmatch('multiarray_1d',structure.type)), return; end
  
  disp([ 'Loading ' structure.filename ' (extracting McCode scan steps)' ]);
  
  % extract all column labels (each word is reversed to ease _I and _ERR search)
  column_labels = flipud(strread(fliplr(structure.variables),'%s'));  % reverse string so that _I and _ERR start words
  
  % get indices that refer to monitor Intensity and Error bars
  monitor_I     = strmatch(fliplr('_I'), column_labels);    % _I
  monitor_ERR   = sort([ strmatch(fliplr('_ERR'), column_labels) strmatch(fliplr('_Err'), column_labels) ]);  % _ERR

  % check consistency
  if (length(monitor_I) ~= length(monitor_ERR))
    warning('McPlot:ScanFormatError', ...
      'File: %s: Found %d Monitor intensity, and %d corresponding error entries', ...
      structure.filename, length(monitor_I), length(monitor_ERR));
  end
  
  % get indices that refer to scan variable parameters (i.e. not monitors)
  scan_labels  = setdiff(1:length(column_labels), [ monitor_I ; monitor_ERR ]);
  xvars        = ~isempty(strmatch(fliplr(strtok(structure.xvars)), column_labels, 'exact'));
  if isempty(xvars)
    warning('McPlot:ScanFormatError', ...
      'File: %s: Can not find scanned variable ''%s'' within columns\n%s', ...
      structure.filename, structure.xvars, structure.variables);
    xvars = transpose(linspace(structure.xlimits(1), structure.xlimits(2), structure.size(1)));
  else
    xvars          = structure.data(:,xvars);
  end
  % generate monitor entries from the scanned data
  for index=1:length(monitor_I)
    this_monitor = fliplr(column_labels{monitor_I(index)});
    this_data = structure;   % initiate single monitor data set from the scan structure
    this_data.type = 'array_1d';
    this_data.size = [ structure.size(1) 4 ];
    this_data.title= [ structure.title ':' this_monitor ];
    this_data.component=this_monitor;
    
    intensity      = structure.data(:,monitor_I(index));
    errors         = structure.data(:,monitor_ERR(index));
    events         = ones(size(intensity));
    this_data.data = intensity;
    this_data.x    = xvars;
    this_data.errors=errors;
    this_data.events=events;
    data = { data{:} this_data };
    disp([ 'Loading ' this_data.filename '#' this_monitor ' (scan steps)' ]);
  end
end % mcplot_split_multiarray

% ==============================================================================

