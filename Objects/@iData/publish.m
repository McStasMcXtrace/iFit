function a = publish(a, filename)
% PUBLISH Publish file containing cells to output file
%   PUBLISH(A) exports object A into the current directory.
%   The supporting output files are stored into an "img" subdirectory.
%
%   PUBLISH(A, FILE) exports object A to given HTML file name FILE.
%
%   PUBLISH(A, DIR) exports A into an HTML document in given directory DIR.
%
% Example: b=publish(iData(peaks)); tf=~isempty(dir(b)); delete(b); tf
% Version: $Date$ $Version$ $Author$
% See also iData, iData/save

if nargin < 2, filename = ''; end

a = saveas(a, filename, 'html');

if nargout == 0
  web(a)
end
