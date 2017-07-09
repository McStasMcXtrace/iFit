function f=csvwrite(a, filename)
% f = csvwrite(s, filename, options) : save iData object into a Comma Separated Values (CSV) file
%
%   @iData/csvwrite function to save data sets as a CSV


  if nargin < 2, filename = []; end
 
  f=saveas(a, filename, 'csv');

