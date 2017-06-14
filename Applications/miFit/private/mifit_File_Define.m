function h = mifit_File_Define(fig)
  % File/Enter expression (model/Data)
  
  % the expression will be evaluated from the main workspace. 
  % All variables therein are usable.
  
  % call an inputdlg in non-modal mode.
  % CallBack: answer must be evaluated, and then sent to mifit for import.
  % answer can be:
  %   filename (~isempty(dir(asnwer))) import data/model
  %   other: sent as is to mifit
  
  Callback = [ 'tmp_cb = getfield(get(gcbf,''UserData''),''Answer''); delete(gcbf);' ...
               'if ~isempty(dir(tmp_cb{1})), mifit(tmp_cb); ' ...
               'else eval([ ''mifit('' tmp_cb{1} '')'' ]); end; ' ...
               'clear tmp_cb;' ];
               
  options.CloseRequestFcn = Callback;
  options.Resize          = 'on';
  options.WindowStyle     = 'normal';
  options.Interpreter     = 'tex';
  
  h = inputdlg('Enter an expression to evaluate (filename, iData, iFunc...)', ...
    'miFit: Enter expression to define new Model/Data set', ...
    2, ...
    { 'conv(lorz, gauss)+strline' }, ...
    options);
  
