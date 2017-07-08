function xmlwrite(a, filename)
% f = xmlwrite(s, filename) : save iFunc object into an XML file
%
%   @iFunc/xmlwrite function to save model as a XML


  if nargin < 2, filename = []; end
 
  f=saveas(a, filename, 'xml');

