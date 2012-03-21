function d = methods(a)
% methods(iData): iData web methods documentation
%
%   @iData/methods: web page documentation
%
%  Type <a href="matlab:doc(iData,'Methods')">doc(iData,'Methods')</a> to access the iFit/Methods Documentation.
%
% Version: $Revision: 1.2 $

% EF 23/10/10 iData impementation

d = [ fileparts(which('iData/version')) filesep '..' filesep 'Docs' filesep 'Methods.html' ];
disp(version(iData))
disp('Opening iData methods documentation from ')
disp(d)
web(d);

