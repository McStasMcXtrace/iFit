function d = doc(a, page)
% doc(iFunc): iFunc web page documentation
%
%   @iFunc/doc: web page documentation
%
%     Open the iFunc documentation.
%     doc(iFunc,page) opens a specific documentation page
%
%     doc(iFunc,'Load')
%     doc(iFunc,'Save')
%     doc(iFunc,'Math')
%     doc(iFunc,'Fit')
%     doc(iFunc,'Plot')
%     doc(iFunc,'Methods')
%
% Version: $Revision: 1.7 $

% EF 23/10/10 iFunc impementation
if nargin ==1, page=''; end
if isempty(page), page='index.html'; end
[p,f,e] = fileparts(page);
if isempty(e), e='.html'; end
page = [ f e ];

url = [ ifitpath filesep 'Docs' filesep page ];
disp(version(iFunc))
disp('Opening iFunc documentation from ')
if length(url)
  disp([ '  <a href="matlab:web ' url '">web ' url '</a>' ]);
else
  disp([ '  ' url ]);
end
web(url);

