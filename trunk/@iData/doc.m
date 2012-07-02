function d = doc(a, page)
% doc(iData): iData web page documentation
%
%   @iData/doc: web page documentation
%
%     Open the iData documentation.
%     doc(iData,page) opens a specific documentation page
%
%     doc(iData,'Load')
%     doc(iData,'Save')
%     doc(iData,'Math')
%     doc(iData,'Fit')
%     doc(iData,'Plot')
%     doc(iData,'Methods')
%
% Version: $Revision: 1.6 $

% EF 23/10/10 iData impementation
if nargin ==1, page=''; end
if isempty(page), page='index.html'; end
[p,f,e] = fileparts(page);
if isempty(e), e='.html'; end
page = [ f e ];

d = [ fileparts(which('iData/version')) filesep '..' filesep 'Docs' filesep page ];
disp(version(iData))
disp('Opening iData documentation from ')

if ~isdeployed
  disp([ '  <a href="matlab:web ' d '">web ' d '</a>' ])
  web(d);
else
  disp(d);
  % attempts to guess how to launch the browser
  if ispc
    system([ 'start /D "' ifitpath '" Docs'])
  elseif ismac
    system([ 'open "' d '"' ]);
  else
    system([ 'gnome-open "' d '"' ]);
  end
end

