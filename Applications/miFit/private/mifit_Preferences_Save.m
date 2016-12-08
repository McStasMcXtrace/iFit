function mifit_Preferences_Save(config)
% File/Preferences: save the Preferences 'config' structure
  filename = fullfile(prefdir, [ 'mifit' '.ini' ]);
  NL = sprintf('\n');
  description = 'miFit interface to iFit';
    str = [ '% miFit configuration script file: ' description NL ...
          '%' NL ...
          '% Matlab ' version ' m-file ' filename NL ...
          '% generated automatically on ' datestr(now) ' with ifit.mccode.org ' 'mifit' NL...
          class2str('config', config) ];
  [fid, message]=fopen(filename,'w+');
  if fid == -1
    warning([ datestr(now) ': Error opening file ' filename ' to save ' description ' configuration.' ]);
    filename = [];
  else
    fprintf(fid, '%s', str);
    fclose(fid);
    mifit_disp([ '[Save_Preferences] Saving Preferences into ' filename ]);
  end
