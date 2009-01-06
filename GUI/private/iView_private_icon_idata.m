function data = iView_private_icon_idata(data, indexes)
% iView_private_icon_idata handles indexes in an iData array
  if ~length(indexes), data=[]; return; end
  valid = find(indexes > 0 & indexes <= length(data));
  indexes = indexes(valid);
  if length(data) > 1, 
    data=data(indexes);
  else
    if length(data) == 0, data=[]; 
    end
  end    
