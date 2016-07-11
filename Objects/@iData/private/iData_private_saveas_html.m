function filename = iData_private_saveas_html(a, filename)

  [Path, name, ext] = fileparts(filename);
  target = fullfile(Path, name);
  titl = char(a);
  titl(titl=='<')='[';
  titl(titl=='>')=']';
  r = report_generator(char(a), target);
  r.open();
  r.section(titl);
  
  % get Model information
  m = []; mp = []; mv = [];
  try
    if isfield(a, 'Model')
      m = get(a, 'Model');
    elseif ~isempty(findfield(a, 'Model'))
      m = get(a, findfield(a, 'Model', 'cache first'));
    end
    
    if isfield(a, 'modelValue')
      mv = get(a, 'modelValue');
    elseif ~isempty(findfield(a, 'modelValue'))
      mv = get(a, findfield(a, 'modelValue', 'cache first'));
    end
    
    if isa(m, 'iFunc')
      % get the parameter values as a struct
      mp = cell2struct(num2cell(m.ParameterValues(:)),strtok(m.Parameters(:)));
    end
  end
  
  % build the HTML content
  r.subsection('Data set')
  % object description
  if ~isempty(a.Title) || ~isempty(title(a))
    data.Title  = strtrim([ a.Title ' ' title(a) ]); end
  if ~isempty(a.Label) || ~isempty(a.DisplayName), 
    data.Label  = strtrim([ a.Label ' ' a.DisplayName ]); end
  data.Source = a.Source;
  data.Date   = [ datestr(a.Date) ', modified ' datestr(a.ModificationDate) ];
  desc = evalc('disp(data)');
  r.add_text([ '<pre> ' desc ' </pre>' ]);
  r.end_section();
  % model description
  if ~isempty(m)
    r.subsection('Model')
    r.add_text(m.Name)
    desc = evalc('disp(mp)');
    r.add_text([ '<pre> ' desc ' </pre>' ]);
    r.end_section();
  end
  
  % plot of the object: special case for 1D which can overlay data and model
  f = [];
  f = figure('Visible','off');
  if ndims(a) == 1 && ~isempty(m) && isa(m, 'iFunc')
    % 1D plot with model
    h=plot(a,'bo',m,'r-','tight');
  elseif ndims(a) == 1
    % simple plot. Model not available.
    h=plot(a,'bo','tight');
  elseif ndims(a) > 1
    if ~isempty(m) && isa(m, 'iFunc')
      h=subplot([a m], [1 2], 'tight');
    else
      h=plot(a, 'tight');
    end
  end
  if ishandle(f)
    set(f, 'Name', [ 'iFit_DataSet_' a.Tag ]);
    r.add_figure(f, char(a), 'centered');
  end
  saveas(f, fullfile(target, 'img', [ 'iFit_DataSet_' a.Tag '.fig' ]), 'fig');
  close(f);
  % export object into a number of usable formats
  builtin('save', fullfile(target, 'img', [ 'iFit_DataSet_' a.Tag '.mat' ]), 'a');
  save(a, fullfile(target, 'img', [ 'iFit_DataSet_' a.Tag '.pdf' ]));
  
  save(a, fullfile(target, 'img', [ 'iFit_DataSet_' a.Tag '.h5' ]), 'mantid');
  if prod(size(a)) < 1e5
    save(a, fullfile(target, 'img', [ 'iFit_DataSet_' a.Tag '.dat' ]), 'dat data');
    save(a, fullfile(target, 'img', [ 'iFit_DataSet_' a.Tag '.svg' ]));
  end
  t = 'Exported to: [ ';
  for ext={'mat','png','dat','svg','pdf','fig','h5'}
    f = [ 'iFit_DataSet_' a.Tag '.' ext{1} ];
    if ~isempty(dir(fullfile(target, 'img', f)))
      t = [ t '<a href="img/' f '"> ' ext{1} ' </a>' ];
    end
  end
  t = [ t ' ]' ];
  r.add_text(t);
  
  % add more information about the Data set
  r.subsection('More about the Data set...')
  desc = evalc('disp(a,''data'',''flat'')');
  desc(desc=='<')='[';
  desc(desc=='>')=']';
  r.add_text([ '<br><pre> ' desc ' </pre>' ]);
  
  % display a 'footer' below the object description
  r.add_text('<hr>');
  r.add_text([ '<b>' datestr(now) '</b> - ' version(iData) '<br>' ])
  
  r.add_text([ '<a href="http://ifit.mccode.org">Powered by iFit ' ...
    '<img src="http://ifit.mccode.org/images/iFit-logo.png" width=35 height=32></a> ' ...
    '<a href="http://www.ill.eu">(c) ILL ' ...
    '<img title="ILL, Grenoble, France www.ill.eu" src="http://ifit.mccode.org/images/ILL-web-jpeg.jpg" alt="ILL, Grenoble, France www.ill.eu" style="width: 21px; height: 20px;"></a><hr>' ]);
  r.add_text('<p><!-- pagebreak --></p> '); % force page break in case we append new stuff
  r.end_section();
  r.close();
  filename = fullfile(target,'index.html');
