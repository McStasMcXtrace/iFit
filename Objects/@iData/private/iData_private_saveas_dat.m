function filename = iData_private_saveas_dat(a, filename)
  NL = sprintf('\n');
  str = [ '# Format: data with text headers' NL ...
            '# URL: http://ifit.mccode.org' NL ...
            '# Creator: iFit/@iData/saveas - ' version(a) NL ...
            '# Title: ' a.Title NL ...
            '# Label: ' a.Label NL ...
            '# DisplayName: ' a.DisplayName NL ...
            '# User: ' a.User NL ...
            '# CreationDate: ' get(a,'Date') NL ...
            '# ModificationDate: ' get(a,'ModificationDate') NL ...
            '# Tag: ' a.Tag NL ];
  n = sprintf('# Type: %dD', ndims(a));
  if ~isvector(a) % histogram
    n = [ n ' data set ' mat2str(size(a)) ];
    str = [ str n NL class2str('', a, 'flat') ];
  else            % event list
    dat = zeros(prod(size(a)), ndims(a)+1);
    lab = '# Data:';
    for index=1:ndims(a)
      dat(:,index) = getaxis(a, index); % store axes first
      this_lab = label(a, index);
      if isempty(this_lab), this_lab=getaxis(a, num2str(index)); end
      lab = [ lab ' ' this_lab ];
    end
    dat(:, ndims(a)+1) = getaxis(a, 0);
    lab = [ lab ' ' 'Signal' ];
    n = [ n ' event list ' mat2str(size(a)) ];
    str = [ str n NL lab NL class2str('', dat, 'flat') ];
  end
  [fid, message]=fopen(filename,'w+');
  if fid == -1
    iData_private_warning(mfilename,[ 'Error opening file ' filename ' to save object ' a.Tag 'in format ' format ]);
    disp(message)
    return
  end
  fprintf(fid, '%s', str);
  fclose(fid);
