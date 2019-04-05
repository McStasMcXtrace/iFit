function c = cell(s)
%  CELL Convert structure array to cell array.
%    C = CELL(S) converts the structure S (with P fields)
%    into a cell array C.
% 
% Example: s = estruct(1:12); iscell(cell(s))
% Version: $Date$ $Author$
% see also estruct, estruct.double, estruct.char

c = struct2cell(s);
