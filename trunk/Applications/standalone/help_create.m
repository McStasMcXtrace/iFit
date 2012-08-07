% create the help strings, to store in .txt files in the Docs/Help directory
function help_create
  to_parse = {'Objects/@iData','Objects/@iFunc','Libraries/Loaders','Libraries/Models','Libraries/Optimizers','Scripts/load','Applications/standalone','Tests'};
  pw = ifitpath;
  disp('Creating help pages for deployed version');
  for index=1:length(to_parse)
    d = to_parse{index};
    disp(d);
    cd (pw);
    cd (d);
    c = dir('*.m');       % get the contents m files
    for fun= 1:length(c); % scan functions
      [p,f,e] = fileparts(c(fun).name);
      h       = help([ d filesep f '.m' ]);
      fid     = fopen([ f '.txt' ],'w+');
      fprintf(fid,'File %s%s%s\n\n', d, filesep, c(fun).name);
      fprintf(fid,'%s\n', h);
      fclose(fid);
    end
  end
