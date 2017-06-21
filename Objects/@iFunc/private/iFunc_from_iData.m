function y = iFunc_from_iData(d)
  % converts an iData object into an iFunc one
  
  % create parameters:
  % Intensity scaling
  % a Scaling and Offset parameter per axis
  %
  % the Expression should interpolate the iData on axes x,y,z,...
  
  y = [];
  
  if numel(d) > 1
    for index=1:numel(d);
      y = [ y feval(mfilename, d(index)) ];
    end
    return
  end
  
  if isempty(d), return; end
  if ischar(d) && ~isempty(dir(d)), d = iData(d); end
  
  y.Description     = sprintf([ 'A model build from the data set\n' ...
    '  ' char(d) '\n\n' ]);
  y.Parameters      = { 'Intensity_scaling' };
  y.ParameterValues = [ 1 ];
  y.Dimension       = ndims(d);
  
  ax = ',x,y,z,t';  % up to 4 dimensions
  axes_when_nan = '';
  axes_homothetic='';
  for index=1:ndims(d)
    lab = label(d, index);
    if isempty(lab), lab = sprintf('Axis%d', index); end
    x = getaxis(d, index); x=x(:);
    y.Parameters{end+1} = [ lab '_offset'  ];
    y.Parameters{end+1} = [ lab '_scaling' ];
    y.ParameterValues(end+1) = 0;
    y.ParameterValues(end+1) = 1;
    y.Description = [ y.Description ...
      sprintf('\n  Axis %i "%s" label is "%s", range [%g:%g]', index, ax(2*index), label(d, index), min(x), max(x)) ];
    axes_when_nan = [ axes_when_nan   ...
      sprintf('%s=getaxis(signal,%i); ', ax(2*index),index) ];
    axes_homothetic=[ axes_homothetic ...
      sprintf('%s=p(%i)+p(%i)*%s; ', ax(2*index), 2*index,2*index+1,ax(2*index)) ];
  end
  
  y.Name          = [ 'iFunc(' char(d) ') model from ' d.Tag ];
  
  y.UserData.Data = d;
  y.Guess = y.ParameterValues;
  
  y.Expression = { 'signal = this.UserData.Data;' };  % must be an iData object
  y.Expression{end+1} = ...
    sprintf([ 'if ~all(isempty(x(:))) && ~all(isnan(x(:)))\n' ...
              '  ' axes_homothetic '\n' ...
              '  signal = interp(signal ' ax(1:(2*y.Dimension)) ');\n' ...
              'else ' axes_when_nan ' end;' ]);
  y.Expression{end+1} = 'signal = double(signal)*p(1);';
  
  y = iFunc(y);
  
  disp([ mfilename ': built model ' y.Name ]);
  
  
