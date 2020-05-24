function f=xmlwrite(a, filename)
% XMLWRITE Write object to an XML Document.
%   F = XMLWRITE(S, FILENAME) saves object S into XML file FILENAME.
%
% Example: f=xmlwrite(iData(peaks)); tf=~isempty(dir(f)); delete(f); tf
% Version: $Date$ $Version$ $Author$
% See also iData, iData/saveas


  if nargin < 2, filename = []; end
 
  f=saveas(a, filename, 'xml');

