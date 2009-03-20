function iView_private_sort(instance, selection)
% iView_private_sort sorts data sets in an iView instance
% selection may be: Size, Date, Title, Label, Tag
  Data = getappdata(instance, 'Data');
	if nargin < 2
		selection = '';
	end
	if isempty(selection)
		items = {'Title (Name)','Size','Date','Label','Tag (unique ID)'};
    selection = listdlg('PromptString', {'Select sorting method to arrange data sets',[ 'in the iView window ' num2str(instance) ]}, ...
      'Name', 'iView: Sort data sets by...', ...
      'SelectionMode','single',...
      'OKString','Sort', 'ListSize', [160 150], ...
      'ListString', items);
    if isempty(selection), return; end
    Data = getappdata(gcf, 'Data');
    selection = items{selection};
  end
  switch selection
  case 'Size'
    list = zeros(1,length(Data));
    for index=1:length(Data)
      list(index) = prod(size(Data(index)));
    end
  case 'Date'
    listd = get(Data,strtok(selection));
    for index=1:length(Data)
      list(index) = datenum(listd(index));
    end
  otherwise
    list = get(Data,strtok(selection));
  end
  [dummy, index] = sortrows(list(:));
  Data = Data(index);
  setappdata(instance, 'Data', Data);
  [hIcon,config,Data]=iView_private_icon(instance, 'check', Data);
  iView_private_documents(instance);
