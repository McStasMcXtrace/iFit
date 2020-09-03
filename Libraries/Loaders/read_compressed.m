function s = read_compressed(filename, option, getpath)
% READ_COMPRESSED Read a compressed file (archive). 
%   Read a ZIP, TAR, GZIP, LZ4, ZSTD, BZIP2, XZ, LZO, ... file.
%   You should install the extractors individually.
%   Recommended:
%   * LZ4 <https://lz4.github.io/lz4/> (2011). Extremely fast.
%   * ZSTD <https://github.com/facebook/zstd> (2015). Very fast.
%   * BROTLI <https://github.com/google/brotli> (2013). Very fast.
%   * LZO <https://www.lzop.org/> (1996-2017). Very fast.
%   * PBZIP2 <https://github.com/ruanhuabin/pbzip2>. Parallelized BZIP2.
%   * PIGZ <https://zlib.net/pigz/> (2007). Parallelized GZIP.
%   * PIXZ <https://github.com/vasi/pixz> (2010). Parallelized XZ.
%
%   Other:
%   * ZIP <https://support.pkware.com/home> (1989). Standard.
%   * TAR <https://www.gnu.org/software/tar/>. Standard.
%   * GZIP <https://www.gnu.org/software/gzip/> (1992). Standard.
%   * BZIP2 <https://www.sourceware.org/bzip2/> (1996). Very compact.
%   * LZMA <https://tukaani.org/lzma/> (1998). Slow, efficient.
%   * XZ <https://tukaani.org/xz/>. Rather slow, very compact.
%   * PXZ <https://jnovy.fedorapeople.org/pxz/>. Parallelized XZ.
%   * 7Z <https://www.7-zip.org/> (1998).
%   * RAR <https://www.rarlab.com/> (1993).
%   * COMPRESS (Z) <https://ncompress.sourceforge.io/> (1985).
%
%   Recommended compressors are LZ4, ZSTD, PIGZ and PBZIP2.
%   ZIP, GZIP and TAR are supported without further installation.
%
%   data=READ_COMPRESSED(file) extract compressed data from given file. When 
%   available, a ramdisk (/dev/shm on Linux) is used as temporary directory.
%
%   data=READ_COMPRESSED(file, tmpdir) extract compressed data from given file, 
%   in given temporary directory. It is recommended to use a RAMdisk for higher
%   efficiency. Temporary files are removed afterwards.
%
%   files=READ_COMPRESSED(..., 'path') only returns the extracted file list. 
%   Data files are NOT imported, and kept decompressed on disk. This about 
%   removing them afterwards, especially when stored in a RAMdisk.

% Not yet done, but possible:
%   READ_COMPRESSED(file, output) compress the data or file into given
%   output filename. The output file extension specifies the archive format.

persistent present

if isempty(present)
  present = check_compressors;
end

s=[];

if nargin == 0 || any(strcmp(filename, {'identify','query','defaults'}))
  % fill in all supported formats
  s = {};
  for index = 1:numel(present)
    s{end+1} = present(index);
  end
  return
end

if nargin < 2, option  = ''; end
if nargin < 3, getpath = ''; end
if isempty(dir(filename)), return; end

if any(strcmp(option,{'path','files','names','list'}))
  getpath = 'path';
  option  = '';
end

% test /run/shm for extraction of temporary files.
if isempty(option)
  for tmp = {'/dev/shm','/run/shm'}
    tf = false;
    if isdir(tmp{1})
      % folder exists. Check free space on disk
      FileObj    = java.io.File(tmp{1});
      free_bytes = FileObj.getFreeSpace;
      delete(FileObj);
      % we need at least 15 times the initial file size (for extraction).
      sz = getfield(dir(filename), 'bytes');
      if sz*15 < free_bytes, tf = true; end
    end
    if tf, option = tmp{1}; break; end
  end
end
% or use system TMP
if isempty(option), option= tempdir; end

if ischar(option) && ~isempty(dir(option)) && ~isdir(option)
  option = fullparts(option); % get path
end

% extract, then call iLoad with resulting files/dir.
if isdir(option)

  % create subdir
  tmpdir = tempname(option);
  [success, message] = mkdir(tmpdir);
  if ~success
    disp([ mfilename ': can not create directory into ' option ' for ' upper(compressor.extension) ]);
    disp(message);
    return; % usually 'permission denied'
  end
  
  try
    % extract
    filenames = decompress(filename, tmpdir, present);
      
    % read content of extracted archive with iLoad
    if ~isempty(filenames)
      if any(strcmp(getpath,{'path','files','names','list'}))
        if numel(filenames) == 1, filenames = filenames{1}; end
        s = filenames;
      else
        s = iLoad(filenames);
      end
    end
  catch ME
    disp([ mfilename ': could not extract ' filename ' into ' option ]);
    disp(getReport(ME));
  end
  
  % remove temporary dir/file
  if isdir(tmpdir) && ~any(strcmp(getpath,{'path','files','names','list'}))
    try
      rmdir(tmpdir, 's');
    end
  end
  
end

% ------------------------------------------------------------------------------
function filenames = decompress(filename, tmpdir, present)
% DECOMPRESS decompress filename locally

  filenames = [];

  % identify which compressor is used (from extension)
  [p,f,e] = fileparts(filename);
  index   = find(strcmp(lower(e(2:end)), lower({ present.extension })));
  if isempty(index), return; end % not an archive
  compressor = present(index);

  % extract it
  if isa(compressor.decompress, 'function_handle') ...
  || strncmp(compressor.decompress, 'matlab:', 7)
    % builtin Matlab extractor
    if strncmp(compressor.decompress, 'matlab:', 7)
      compressor.decompress = compressor.decompress(8:end);
    end
    filenames = feval(compressor.decompress, filename, tmpdir);
  else
    % system call
    % copy initial file in temporary directory
    copyfile(filename, fullfile(tmpdir, [f e]));
    
    % extract
    system([compressor.decompress ' ' fullfile(tmpdir, [f e]) ]);
    
    % remove initial archive from temp dir
    if ~isempty(dir(fullfile(tmpdir, [f e])))
      delete(fullfile(tmpdir, [f e]));
    end
    
    filenames = tmpdir; % we shall extract all files in temp dir (from archive)

  end
  
% ------------------------------------------------------------------------------
function present = check_compressors(options)
% CHECK_COMPRESSORS check if (de)compressors are present
%
% Supported:
%   ZIP, GZIP, TAR                    (builtin)
%   LZ4, COMPRESS/Z, ZSTD, LZO, BZIP2 (via system)

  present = [];

  % required to avoid Matlab to use its own libraries
  if ismac,      precmd = 'DYLD_LIBRARY_PATH= ; DISPLAY= ; ';
  elseif isunix, precmd = 'LD_LIBRARY_PATH= ;  DISPLAY= ; ';
  else           precmd = ''; end

  %         { EXT,  CHECK,          cmd UNCOMP,      cmd COMP}
  % last extractors should be the default (and less efficient)
  tocheck = {'lz4', 'lz4 --version',    'lz4 -d',        'lz4 -1';
             'zst', 'zstd --version',   'zstd -d',       'zstd -1',;
             'Z',   'compress -V',      'compress -d',    'compress';
             'lzo', 'lzop --version',   'lzop -d',       'lzop -1';
             'gz',  'pigz --version',   'pigz -d',       'pigz -1';
             'bz2', 'pbzip2 --version', 'pbzip2 -d',     'pbzip2 -1';
             'bz2', 'bzip2 --version',  'bzip2 -d',      'bzip2 -1';
             'br',  'brotli --version', 'brotli -d',     'brotli -1';
             'xz',  'pixz -version',    'pixz -d',       'pixz -1';
             'xz',  'pxz --version',    'pxz -d',        'pxz -1';
             'xz',  'xz --version',     'xz -d',         'xz -1';
             'rar', 'rar -version',     'rar x',         'rar a';
             '7z',  '7z i',             '7z x',          '7z a';
             'lzma','lzma --version',   'lzma -d',       'lzma -1';
             'zip', 'matlab:unzip',     'matlab:unzip',  'matlab:zip';
             'gz',  'matlab:gunzip',    'matlab:gunzip', 'matlab:gzip';
             'tar', 'matlab:untar',     'matlab:untar',  'matlab:tar' };

  for totest = tocheck'
    % look for executable and test with various extensions
    ok = false;
    if isa(totest{2}, 'function_handle') ...
    || strncmp(totest{2}, 'matlab:', 7)
      ok = true;
    else
      [status, result] = system([ precmd totest{2} ]); % usually 127 indicates 'command not found'
    end
    if (isunix && any(status == [0 2])) ...
    || (ispc &&  isempty(strfind(result, [ '''' strtok(totest{1}) '''' ])))
      ok = true;
    end
    if ok && (~isstruct(present) || ~any(strcmp(totest{1}, { present.extension })))
      p.name        = [ 'Compressed archive ' totest{1} ];
      p.extension   = totest{1};
      p.method      = mfilename;
      p.decompress  = totest{3};
      p.compress    = totest{4};
      if isstruct(present)
        present(end+1) = p;
      else
        present = p;
      end
    end
  end
