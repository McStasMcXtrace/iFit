function d = doc(a)
% doc(iData): iData web page documentation
%
%   @iData/doc: web page documentation
%
% Version: $Revision: 1.1 $

% EF 23/10/10 iData impementation

d = [ fileparts(which('iData/version')) filesep '..' filesep 'Docs' filesep 'index.html' ];
disp('Opening iData documentation from ')
disp(d)
web(d);

