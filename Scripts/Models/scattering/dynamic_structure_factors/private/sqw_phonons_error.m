function sqw_phonons_error(message, options)

if options.gui && ishandle(options.gui)
  delete(options.gui);
  errordlg(message, [ 'iFit: ' mfilename ' ' options.configuration ' FAILED' ]);
end
if ~isdeployed && usejava('jvm') && usejava('desktop')
  disp([ '<a href="matlab:doc(''' mfilename ''')">help ' mfilename '</a> (click here to get help)' ])
end
sqw_phonons_htmlreport(fullfile(options.target, 'index.html'), 'error', options, message);
error(message);
