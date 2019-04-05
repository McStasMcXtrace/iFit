function out = load(a, filename)
% d = load(s, file): iFunc Model loader
%
%   @iFunc/load: imports a Model from a file into Matlab/iFunc object
%   The input argument 'file' should be a file name or empty (then popup a file selector).
%
%   Default supported formats include: MAT, M, YAML, JSON, XML format holding a model (iFunc/save)
%   The list of supported formats to create iFunc objects is obtained with:
%     load(iFunc,'formats')
%
% input:  s: object (iFunc)
%         file: file name(s) to import (char/cellstr)
%
% output: d: single object or array (iFunc)
% ex:     save(gauss, 'gauss.yaml'); load(iFunc, 'gauss.yaml');
%
% Version: $Date$ $Version$ $Author$
% See also: iLoad, save, iFunc/save, iData/save, Loaders
  out = [];
  if nargin < 2
    filename = '';
  end
  
  % supported format list
  filterspec = {...
      '*.*','All files (*.*)'; ...
      '*.json', 'JSON JavaScript Object Notation (*.json)'; ...
      '*.m',   'Matlab script/function (*.m)'; ...
      '*.mat', 'Matlab binary file (*.mat)'; ...
      '*.mat', 'SpinW object stored into a Matlab binary file (*.mat)'; ...
      '*.xml','XML file (*.xml)'; ...
      '*.yaml;*.yml','YAML interchange format (*.yaml)' };

  if strcmp(filename, 'formats')
    fprintf(1, '       EXT  DESCRIPTION [%s(iFunc)]\n', mfilename);
    fprintf(1, '-----------------------------------------------------------------\n'); 
    for index=1:size(filterspec,1)
      ext = upper(filterspec{index,1});
      ext = strrep(ext,'.','');
      ext = strrep(ext,'*','');
      fprintf(1,'%10s  %s \n', ext, filterspec{index,2});
    end
    out = filterspec;
    return
  end
  if isempty(filename)
    % pop-up  GUI to select a Model file
    [filename, pathname, filterindex] = uigetfile( ...
       filterspec, ...
        ['Load iFunc/Model from ...']);
    if ~isempty(filename) & filename ~= 0
      ext = filterspec{filterindex,1};
      if iscell(ext) && ischar(ext{1}), ext=ext{1}; end
      % check if extension was given
      [f,p,e] = fileparts(filename);
      if isempty(e), 
        filename=[ filename ext(2:end) ];
        format=ext(3:end);
      else
        format=e(2:end);
      end
      filename = strcat(pathname, filesep, filename);
    else
      return
    end
  end

% import with the selected filename
  out = iFunc(filename);
