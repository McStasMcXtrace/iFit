function c = cell(s)
%  CELL Convert structure array to cell array.
%    C = CELL(S) converts the structure S (with P fields)
%    into a cell array C.
% 
% Example: s = iData(1:12); iscell(cell(s))
% Version: $Date$ $Version$ $Author$
% see also iData, iData.double, iData.char

c = struct2cell(s);
