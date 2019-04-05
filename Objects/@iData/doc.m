function url = doc(a, token)
% doc(iData): iData web page documentation
%
%   @iData/doc: web page documentation
%
%     Open the iData documentation.
%     doc(iData, page) opens a specific documentation page
%
%   You may as well search for a token with:
%     doc(iData, token)
%
%     doc(iData,'Load')
%     doc(iData,'Save')
%     doc(iData,'Math')
%     doc(iData,'Fit')
%     doc(iData,'Plot')
%     doc(iData,'Methods')
%     doc(iData,'Sqw2D')
%
% Version: $Date$ $Version$ $Author$

% EF 23/10/10 iData implementation
if nargin ==1, token=''; end
if isempty(token), token='index.html'; end
[p,f,e] = fileparts(token);
if isempty(e), e='.html'; end
page = [ f e ];

url      = fullfile(ifitpath,'Docs',page);
url_path = strtok(url, '#?');

if isempty(dir(url_path)) % page does not exist ?
  % search for token in whole documentation
  [a,b]=grep('-i','-n','-r', token, ...
    { fullfile(ifitpath,'Docs'),   fullfile(ifitpath,'Scripts') ...
      fullfile(ifitpath,'Objects') fullfile(ifitpath,'Libraries') ...
      fullfile(ifitpath,'Applications')});
  
  % a contains the file names
  % b.findex lists occurrences per file = b.lcount
  % we sort the matching file entries (a) as a function of the match occurrences
  [~,index] = sort(b.lcount,1,'descend');
  url = a(index);         % files sorted by decreasing occurrences
  occ = b.lcount(index);  % occurrences, decreasing
  
  % now sort results in: html, txt, m, then the rest
  i_html=[]; i_txt=[]; i_m=[]; i_other=[];
  for index=1:numel(url)
    [p,f,e] = fileparts(url{index});
    switch lower(e)
    case {'.html','.htm'}
      i_html(end+1) = index;
    case {'.txt'}
      i_txt(end+1) = index;
    case {'.m'}
      i_m(end+1) = index;
    otherwise
      i_other(end+1) = index;
    end
  end
  index= [ i_html i_txt i_m i_other ];
  url = url(index);
  occ = occ(index);
  
  if isempty(a)
      disp([ mfilename ': token ' token ' not found in Help.' ]);
      return; 
  end
  if numel(url) > 1
    % display a listdlg([ url occ ]) to select which file to open
    matches   = strcat('(',cellstr(num2str(occ)),') ',url);
    selection = listdlg('ListString', matches, 'ListSize', [ 300 160 ], ...
      'Name', [ mfilename ': iFit Help pages mentioning ' token ], ...
      'SelectionMode', 'single', ...
      'PromptString', { [ 'Here are the iFit entries for "'  token '".']; ...
      'Showing occurrences and help resource, sorted as HTML,text,m and other files.'; ...
      'Please choose one of them to display.' });
    if isempty(selection),    return; end % cancel
  else selection = 1;
  end
  if iscell(url) url = url{selection}; end
end

disp(version(iData))
disp('Opening iData documentation from ')
if length(url) && ~isdeployed && usejava('jvm')
  disp([ '  <a href="matlab:web ' url '">web ' url '</a>' ]);
else
  disp([ '  ' url ]);
end
[p,f,e] = fileparts(url);
switch lower(e);
case {'.txt','.m'}
  if isdeployed
    TextEdit(url);
  else
    edit(url);
  end
otherwise
  web(url);
end

