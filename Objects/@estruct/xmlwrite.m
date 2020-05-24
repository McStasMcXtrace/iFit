function f=xmlwrite(a, filename)
% XMLWRITE Write object to an XML Document.
%   F = XMLWRITE(S, FILENAME) saves object S into XML file FILENAME.
%
% Example: f=xmlwrite(estruct(peaks)); ~isempty(dir(f)); delete(f);
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/saveas


  if nargin < 2, filename = []; end
 
  f=saveas(a, filename, 'xml');

