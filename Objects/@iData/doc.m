function url = doc(a, page)
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
% Version: $Date$

% EF 23/10/10 iData impementation
if nargin ==1, page=''; end
if isempty(page), page='index.html'; end
[p,f,e] = fileparts(page);
if isempty(e), e='.html'; end
page = [ f e ];

url = fullfile(ifitpath,'Docs',page);
url_path = strtok(url, '#?');
if ~isempty(dir(url_path)) % page exists ?
  disp(version(iData))
  disp('Opening iData documentation from ')
  if length(url) && ~isdeployed && usejava('jvm')
    disp([ '  <a href="matlab:web ' url '">web ' url '</a>' ]);
  else
    disp([ '  ' url ]);
  end
  web(url);
else
  % search for token in whole documentation
  [a,b]=grep('-i','-n','-r', page, fullfile(ifitpath,'Docs'));
  % a contains the file names
  % b.findex lists occurencies per file = b.lcount
  
  % we sort the matching file entries (a) as a function of the match occurencies
  [~,index] = sort(b.lcount,1,'descend');
  url = a(index);         % files sorted by decreasing occurencies
  occ = b.lcount(index);  % occurencies, decreasing
  
  if numel(url) > 1
    % display a listdlg([ url occ ]) to select which file to open
    matches   = strcat('(',cellstr(num2str(occ)),') ',url);
    selection = listdlg('ListString', matches, 'ListSize', [ 300 160 ], ...
      'Name', [ mfilename ': iFit Help pages mentioning ' page ], ...
      'SelectionMode', 'single', ...
      'PromptString', { [ 'Here are the iFit entries for ' page ]; ...
      'Showing number of matches and help resource.'; ...
      'Please choose one of them to display.' });
    if isempty(selection),    return; end % cancel
  else selection = 1;
  end
  if iscell(url), url = url{selection}; end
  web(url);
end

