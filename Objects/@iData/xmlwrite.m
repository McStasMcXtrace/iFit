function f=xmlwrite(a, filename)
% f = xmlwrite(s, filename) : save iData object into an XML file
%
%   @iData/xmlwrite function to save data sets as a XML


  if nargin < 2, filename = []; end
 
  f=saveas(a, filename, 'xml');

