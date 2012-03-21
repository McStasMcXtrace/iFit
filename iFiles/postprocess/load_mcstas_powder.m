function a=load_mcstas_powder(a)
% function a=load_mcstas_powder(a)
%
% Returns an iData style dataset from a McStas Powder file (LAZ/LAU)
% such files can be obtained from Crystallographica and ICSD <icsd.ill.fr>
%
% Version: $Revision: 1.2 $
% See also: iData/load, iLoad, save, iData/saveas

if ~isa(a,'iData')
  a = load(iData,a,'McStas powder');
  return
end

% handle input iData arrays
if length(a(:)) > 1
  for index=1:length(a(:))
    a(index) = feval(mfilename, a(index));
  end
  return
end

a=iData(a);
if isempty(findstr(a,'Lazy')) && isempty(findstr(a,'Crystallographica'))
  warning([ mfilename ': The loaded data set ' a.Tag ' "' a.Title '" is not an Powder text data format.' ]);
  return
end

% set HKL axes and Intensity signal. No Error.
data_definition = getalias(a, 'Signal');
columns_header = getfield(a.Headers, fliplr(strtok(fliplr(data_definition),'.')));
this = findstr(a,'VALUE');
if ~isempty(this), this=this{1}; end
% the header line may be split as it contains numerics. Prepend Headers.VALUE.
columns_header = [ this ' ' columns_header ];
% the Lazy format has a column named 'D VALUE': remove the space so that columns are not shifted
columns_header = strrep(columns_header, 'D VALUE','D_VALUE');
columns = strread(columns_header,'%s','delimiter',' ;#');
columns = columns(~cellfun('isempty', columns));
for index=1:length(columns)
  % clean the column name so that it looks like a variable name
  columns{index} = strrep(columns{index}, '.','');
  columns{index} = strrep(columns{index}, '-','');
  columns{index} = strrep(columns{index}, '/','');
  columns{index} = strrep(columns{index}, '*','');
  columns{index} = strrep(columns{index}, '(','');
  columns{index} = strrep(columns{index}, ')','');
  columns{index} = genvarname(columns{index});
  if ~isfield(a, columns{index})
    setalias(a, columns{index}, [ data_definition '(:,' num2str(index) ')' ]);
    disp([ columns{index} '=' data_definition '(:,' num2str(index) ')' ]);
    this_axis = find(strcmpi(columns{index},{'h','k','l'}),1);
    if ~isempty(this_axis)
      this_axis = this_axis(1);
      setaxis(a, this_axis, columns{index}); 
    end
    this_axis=[];
    if   ~isempty(strfind(columns{index}, 'FHKL')) ...
      || ~isempty(strfind(columns{index}, 'Fsquared'))
      setaxis(a, 'Signal', columns{index});
    end
  end
end


setalias(a,'Error',0);
a = transpose(a);

