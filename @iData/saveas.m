function filename = saveas(a, varargin)
% f = saveas(s, filename, format) : save iData object into various data formats
%
%   @iData/saveas function to save data sets
%   This function save the content of iData objects. 
%
% input:  s: object or array (iData)
%         filename: name of file to save to. Extension, if missing, is appended (char)
%                   If the filename already exists, the file is overwritten.
%         format: data format to use (char), or determined from file name extension
%           'm'   save as a flat Matlab .m file (a function which returns an iData object or structure)
%           'mat' save as a '.mat' binary file (same as 'save')
%           'hdf' save as an HDF5 data set
%           'cdf' save as NetCDF 
%           'gif','bmp' save as an image (no axes, only for 2D data sets)
%           'png','tiff','jpeg','ps','pdf','ill','eps' save as an image (with axes)
%           'xls' save as an Excel sheet (requires Excel to be installed)
%
% output: f: filename used to save file(s) (char)
% ex:     b=saveas(a, 'file', 'm');
%
% See also iData, iData/load, save

if length(a) > 1
  if length(varargin) >= 1, filename_base = varargin{1}; 
  else filename_base = ''; end
  filename = cell(size(a));
  for index=1:length(a(:))
    filename{index} = saveas(a(index), varargin{:});
    if isempty(filename_base), filename_base = filename{index}; end
    if length(a(:)) > 1
      [path, name, ext] = fileparts(filename_base);
      varargin{1} = [ path name '_' num2str(index) ext ];
    end
  end
  return
end

if nargin < 2, filename = ''; else filename = varargin{1}; end
if isempty(filename), filename = a.Tag; end
if nargin < 3, format=''; else format = varargin{2}; end

% handle extensions
[path, name, ext] = fileparts(filename);
if isempty(ext) & ~isempty(format), 
  ext = [ '.' format ]; 
  filename = [ filename ext ];
elseif isempty(format) & ~isempty(ext)
  format = ext(2:end);
elseif isempty(format) & isempty(ext) 
  format='m'; filename = [ filename '.m' ];
end

switch format
case 'jpg'
  format='jpeg';
case 'eps'
  format='epsc';
case 'ps'
  format='psc';
case 'nc'
  format='cdf';
end

switch format
case 'm'
  str = [ 'function this=' name sprintf('\n') ...
          '% Original data: ' class(a) ' ' inputname(1) ' ' a.Tag sprintf('\n') ...
          '%' sprintf('\n') ...
          '% Matlab ' version ' m-file ' filename ' saved on ' datestr(now) ' with iData/saveas' sprintf('\n') ...
          '% To use import data, type ''' name ''' at the matlab prompt.' sprintf('\n') ...
          '% You will obtain an iData object (if you have iData installed) or a structure.' sprintf('\n') ...
          '%' sprintf('\n') ...
          iData_saveas_single('this', a) ];
  [fid, message]=fopen(filename,'w+');
  if fid == -1
    iData_private_warning(mfilename,[ 'Error opening file ' filename ' to save object ' a.Tag ]);
    return
  end
  fprintf(fid, '%s', str);
  fclose(fid);
case 'mat'
  save(filename, a);
case {'hdf','cdf'}
  [fields, types, dims] = findfield(a);
  towrite={};
  for index=1:length(fields(:)) % get all field names
    val=get(a, fields{index});
    if iscellstr(val), 
      val=val(:);
      val(:, 2)={ ';' }; 
      val=val'; 
      val=[ val{1:(end-1)} ];
    end
    if ~isnumeric(val) & ~ischar(val), continue; end
    if strcmp(format,'cdf')
      if isempty(towrite)
        towrite={ fields{index}, val };
      else
        towrite={ towrite{1:end}, fields{index}, val };
      end
    else
      if isempty(towrite)
        hdf5write(filename, [ filesep 'iData' filesep fields{index} ], val);
      else
        hdf5write(filename, [ filesep 'iData' filesep fields{index} ], val, 'WriteMode', 'append');
      end
      towrite=1;
    end
  end
  if strcmp(format,'cdf')
    cdfwrite(filename, towrite);
  end
case 'xls'
  xlswrite(filename, double(a), a.Title);
case 'csv'
  csvwrite(filename, double(a));
case {'gif','bmp','pbm','pcx','pgm','pnm','ppm','ras','xwd'}
  if ndims(a) == 2
    a=double(a);
    b=(a-min(a(:)))/(max(a(:))-min(a(:)))*64;
    imwrite(b, jet(64), filename, format);
  else
    iData_private_warning(mfilename,[ 'Can not save object ' a.Tag ' (non 2D) into image format ' format ]);
    filename='';
  end
case 'epsc'
  f=figure('visible','off');
  plot(a,'view2 axis tight');
  print(f, [ '-depsc -tiff' ], filename);
  close(f);
case {'png','tiff','jpeg','psc','pdf','ill'}
  f=figure('visible','off');
  plot(a,'view2 axis tight');
  print(f, [ '-d' format ], filename);
  close(f);
end

% ============================================================================
% private function
function str=iData_saveas_single(this, data)
% create a string [ 'this = data;' ]

NL = sprintf('\n');
if ischar(data)
  str = [ this ' = ''' iData_saveas_validstr(data) ''';' NL ];
elseif isa(data, 'iData')
  str = [ '% ' this ' (' class(data) ') size ' num2str(size(data)) NL ];
  str = [ str iData_saveas_single(this, struct(data)) ];
  str = [ str 'if ~exist(''iData''), return; end' NL ];
  str = [ str this '_s=' this '; ' this ' = rmfield(' this ',' this '.Alias); ' this ' = iData(' this ');' NL ...
         'setalias(' this ', ' this '_s.Alias.Names, ' this '_s.Alias.Values, ' this '_s.Alias.Labels);' NL ... 
         'setaxis('  this ', mat2str(1:length(' this '_s.Alias.Axis)), ' this '_s.Alias.Axis);' NL ...
         '% end of iData ' this NL ];
elseif isnumeric(data) | islogical(data)
  str = [ '% ' this ' numeric (' class(data) ') size ' num2str(size(data)) NL ...
          this ' = ' mat2str(data(:)) ';' NL this ' = reshape(' this ', ' num2str(size(data)) ');' NL ];
elseif isstruct(data)
  f = fieldnames(data);
  str = [ '% ' this ' (' class(data) ') length ' num2str(length(f)) NL ];
  for index=1:length(f)
    str = [ str iData_saveas_single([ this '.' f{index} ], getfield(data, f{index})) ];
  end
  str = [ str '% end of struct ' this NL ];
elseif iscellstr(data)
  str = [ '% ' this ' (' class(data) 'str) size ' mat2str(size(data)) NL ...
          this ' = { ...' NL ];
  for index=1:length(data(:))
    str = [ str '  ''' iData_saveas_validstr(data{index}) '''' ];
    if index < length(data(:)), str = [ str ', ' ]; end
    str = [ str ' ...' NL ];
  end
  str = [ str '}; ' NL ];
  str = [ str this ' = reshape(' this ',' mat2str(size(data)) ');' NL '% end of cellstr ' this NL ];
elseif iscell(data)
  str = [ '% ' this class(data) ' size ' mat2str(size(data)) NL ...
          this ' = cell(' mat2str(size(data)) ');' NL ];
  for index=1:length(data(:))
    str = [ str iData_saveas_single([ this '{' num2str(index) '}' ], data{index}) NL ];
  end
  str = [ str this ' = reshape(' this ',' mat2str(size(data)) ');' NL '% end of cell ' this NL ];
elseif isa(data, 'function_handle')
  iData_private_warning(mfilename,  'can not save function handles. Skipping.');
else
  iData_private_warning(mfilename,[ 'can not save ' class(data) '. Skipping.' ]);
  % other object
end

function str=iData_saveas_validstr(str)
index = find(str < 32 | str > 127);
str(index) = ' ';

