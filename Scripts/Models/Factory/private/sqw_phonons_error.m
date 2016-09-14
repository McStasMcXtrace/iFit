function sqw_phonons_error(message, options)

if ~isempty(options.gui) && ~isnan(options.gui) && ishandle(options.gui)
  delete(options.gui);
  errordlg(message, [ 'iFit: ' mfilename ' ' options.configuration ' FAILED' ]);
end
if ~isdeployed && usejava('jvm') && usejava('desktop')
  disp([ '<a href="matlab:doc(''sqw_phonons'')">help sqw_phonons</a> (click here to get help)' ])
end
sqw_phonons_htmlreport('', 'error', options, message);
disp(message);
