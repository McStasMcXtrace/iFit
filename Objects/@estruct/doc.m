function url = doc(a, token)
% DOC web page documentation
%
%   DOC(estruct) opens the estruct documentation.
%
%   DOC(estruct, page) opens a specific documentation page.
%
%   DOC(estruct, token) searches for a token.
%
% Example: doc(estruct,'Load'); 1
%
% Version: $Date$ $Version$ $Author$
% See also: estruct.help estruct.edit estruct.inputdlg

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

disp(version(estruct))
disp('Opening estruct documentation from ')
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

