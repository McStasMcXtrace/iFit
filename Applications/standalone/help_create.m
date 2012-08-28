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
  
  return
  
  % now we create launchers for Linux (OpenDesktop .desktop files) and Windows (.bat files)
  
  % categories for launchers ---------------------------------------------------
  
  % Model list (predefined iFunc)
  d = dir([ fileparts(which('gauss')) ]);
  criteria = []; 
  for index=1:length(d)
    this = d(index);
    try
      [dummy, method] = fileparts(this.name);
      options = feval(method,'identify');
      if isa(options, 'iFunc')
        launcher_write(method);
      end
    end
  end % for
  
  % Mathematical operators (iFunc), save 'ans'
  % abs del2 floor sparse transpose  acos conj full sqrt uminus  acosh real  asin exp ndims round xcorr  asinh imag norm  atan cos isempty not sign tan  atanh cosh fliplr log sin tanh  ceil ctranspose flipud log10 sinh minus conv convn xcorr fits
  % Commands (iFunc)
  % edit plot char copyobj doc feval  get subplot 
  
  % Mathematical operators (iData), save 'ans'
  % abs acos asin atan cos sin tan cosh, sinh, tanh acosh, asinh, atan exp log log10 sqrt ctranspose
  % transpose permute floor ceil round sign uminus imag real fft ifft del2 gradient diff sum prod trapz cumsum 
  % plus minus times mtimes rdivide combine power lt le gt ge ne eq conv xcorr interp 
  % mean max mean median std peaks camproj cumtrapz norm cat dog linspace logspace intersect union hist
   
  % Commands (iData)
  % edit plot char copyobj doc feval  get subplot surf mesh load contour surfc surfl plot3 scatter3 waterfall
  % image caxis colormap slice
  
