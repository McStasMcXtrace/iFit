function [filename,format] = save(a, varargin)
% SAVE Save object in desired output format.
%   F =SAVE(S) saves object into a MAT file name determined from the object Tag.
%   The generated filename F is returned.
%   SAVE is an alias to SAVEAS.
%
%   F = SAVE(S,'FILE.EXT') saves the object S into file FILE with format
%   corresponding with extension EXT.
%
%   F = SAVE(S,'FILE','FORMAT') saves the object S into file FILE with format
%   FORMAT. FORMAT can be specified as an extension '.EXT' or 'EXT'.
%
%   SAVE(iData,'formats') prints a list of supported export formats.
%
%  Type <a href="matlab:doc(iData,'Save')">doc(iData,'Save')</a> to access the iFit/Save Documentation.
%
% Example: f=save(iData(peaks)); tf=~isempty(dir(f)); delete(f); tf
% Version: $Date$ $Version$ $Author$
% See also iData, iData/saveas, iData/load, save, iData/plot

[filename,format] = saveas(a, varargin{:});
