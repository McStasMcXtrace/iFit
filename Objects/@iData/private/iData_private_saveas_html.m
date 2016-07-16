function filename = iData_private_saveas_html(a, filename)

  [Path, name, ext] = fileparts(filename);
  target = fullfile(Path, name);
  titl = char(a);
  titl(titl=='<')='[';
  titl(titl=='>')=']';
  % Open and write the HTML header
  filename = fullfile(target,'index.html');
  if ~isdir(target), mkdir(target); end
  if ~isdir(fullfile(target,'img')), mkdir(fullfile(target,'img')); end
  fid = fopen(filename, 'a+');
    
  if fid == -1, filename = []; return; end
  fprintf(fid, '<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n');
  fprintf(fid, '<html>\n<head>\n<title>%s</title>\n<\head>\n', ...
      titl);
  fprintf(fid, '<body><div style="text-align: center;">\n');
  fprintf(fid, '<a href="http://ifit.mccode.org"><img title="ifit.mccode.org" src="http://ifit.mccode.org/images/iFit-logo.png" align="middle" height=100></a>\n');
  fprintf(fid, '<h1>%s</h1></div>\n', titl);
  
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
  fprintf(fid,'<h2>Data set</h2>\n');
  % object description
  if ~isempty(a.Title) || ~isempty(title(a))
    data.Title  = strtrim([ a.Title ' ' title(a) ]); end
  if ~isempty(a.Label) || ~isempty(a.DisplayName), 
    data.Label  = strtrim([ a.Label ' ' a.DisplayName ]); end
  data.Source = a.Source;
  data.Date   = [ datestr(a.Date) ', modified ' datestr(a.ModificationDate) ];
  desc = evalc('disp(data)');
  fprintf(fid,[ '<pre> ' desc ' </pre>\n' ]);

  % model description
  if ~isempty(m)
    fprintf(fid,'<h2>Model</h2>\n');
    fprintf(fid, m.Name);
    desc = evalc('disp(mp)');
    fprintf(fid, [ '<pre> ' desc ' </pre>\n' ]);
  end
  
  % plot of the object: special case for 1D which can overlay data and model
  f = figure('Visible','off', 'Name', [ 'iFit_DataSet_' a.Tag ]);
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
  saveas(f, fullfile(target, 'img', [ 'iFit_DataSet_' a.Tag '.png' ]), 'png');
  fprintf(fid, '<img src="%s" align="middle"><br>\n', ...
    fullfile('img',[ 'iFit_DataSet_' a.Tag '.png' ]));
  saveas(f, fullfile(target, 'img', [ 'iFit_DataSet_' a.Tag '.pdf' ]), 'pdf');
  close(f);
  % export object into a number of usable formats
  builtin('save', fullfile(target, 'img', [ 'iFit_DataSet_' a.Tag '.mat' ]), 'a');

  save(a, fullfile(target, 'img', [ 'iFit_DataSet_' a.Tag '.h5' ]), 'mantid');
  if prod(size(a)) < 1e5
    save(a, fullfile(target, 'img', [ 'iFit_DataSet_' a.Tag '.dat' ]), 'dat data');
    save(a, fullfile(target, 'img', [ 'iFit_DataSet_' a.Tag '.svg' ]));
  end
  t = 'Exported to: [ ';
  for ext={'mat','png','dat','svg','pdf','fig','h5'}
    f = [ 'iFit_DataSet_' a.Tag '.' ext{1} ];
    if ~isempty(dir(fullfile(target, 'img', f)))
      t = [ t '<a href="img/' f '">' ext{1} '</a> ' ];
    end
  end
  t = [ t ' ]' ];
  fprintf(fid, '%s\n', t);
  
  % add more information about the Data set
  fprintf(fid,'<h2>More about the Data set...</h2\n');
  desc = evalc('disp(a,''data'',''flat'')');
  desc(desc=='<')='[';
  desc(desc=='>')=']';
  fprintf(fid,[ '<br><pre> ' desc ' </pre>\n' ]);
  
  % display a 'footer' below the object description
  fprintf(fid,'<hr>\n');
  fprintf(fid,[ '<b>' datestr(now) '</b> - ' version(iData) '<br>\n' ]);
  
  fprintf(fid,[ '<a href="http://ifit.mccode.org">Powered by iFit ' ...
    '<img src="http://ifit.mccode.org/images/iFit-logo.png" width=35 height=32></a> \n' ...
    '<a href="http://www.ill.eu">(c) ILL ' ...
    '<img title="ILL, Grenoble, France www.ill.eu" src="http://ifit.mccode.org/images/ILL-web-jpeg.jpg" alt="ILL, Grenoble, France www.ill.eu" style="width: 21px; height: 20px;"></a><hr>\n' ]);
  fprintf(fid,'<p><!-- pagebreak --></p>\n'); % force page break in case we append new stuff
  fclose(fid);

