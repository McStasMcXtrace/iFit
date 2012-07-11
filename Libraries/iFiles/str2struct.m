function s=str2struct(string)
% s=str2struct(string) Create a structure from string lines
%   This function creates a structure from string containing <name> <value> pairs
%   Structure member must be separated by ';' or end-of-lines.
%   The member assignation can be specified with spaces, '=' and ':'
%
% input arguments:
%   string: string from which the structure should be extracted
%
% output variables:
%   s: structure which contains the named values
%
% example: str=str2struct('Temperature: 200; RV=3; comment something nice');
%          
% See also: mat2str, num2str, eval, sprintf, class2str
%
% Part of: iFiles utilities (ILL library)
% Author:  E. Farhi <farhi@ill.fr>. $Revision: 1.5 $

s={};
if nargin ==0, return; end
if isempty(string), return; end
if ~ischar(string) && ~iscell(string), return; end

% transform the string into a cellstr
string = cellstr(string);

cellstring = {};

% split the string into seperate lines if they contain <EOL> characters
for index=1:numel(string)
  this = string{index};
  split = textscan(this,'%s','Delimiter',sprintf('\n\r\f;'));
  for j=1:numel(split)
    this_split=split{j};
    cellstring = { cellstring{:} this_split{:} };
  end
end

% interpret the line as <name> <separator> <value> <comment>
for index=1:numel(cellstring)
  this = cellstring{index};
  [name, line] = strtok(this, sprintf('=: \t'));
  if isempty(name), continue; end
  if name(1)=='#' || name(1)=='%' || strncmp(name, '//', 2) || name(1) == '!'
    continue; % skip comment lines
  end
  nextline = min(find(isstrprop(line, 'alphanum')));
  startline=line(1:nextline);
  nextline=max(find(startline == '=' | startline == ' ' | startline == ':'));
  if nextline >= 1, nextline=nextline+1; else continue; end
  line = line(nextline:end);
  [value, count, errmsg, nextindex] = sscanf(line, '%f');
  comment = strtrim(line(nextindex:end)); comment(~isstrprop(comment,'print')) = ' ';
  name = strrep(name, '.', '_');
  name = strrep(name, '-', '_');
  name = genvarname(name);
  if isempty(value), 
      value=comment;
      tmp  =str2num(value);
      if ~isempty(tmp) && isnumeric(tmp)
          value=tmp;
      end
  end
  if ~isempty(value), s.(name) = value; end
end


