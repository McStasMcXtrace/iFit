function mifit_Data_Select_ROI(varargin)
% Data/Select_ROI: select a ROI, then create a new selected data set, and its mask as [0-1]

  persistent first_use
  
  if isempty(first_use)
    first_use = false;
    doc(iData, 'Plot.html#mozTocId313946');
  end
  
  D=mifit_List_Data_pull;
  
  mifit_disp([ '[Data select ROI] Selecting ROI on ' char(D{1}) ]);
  set(mifit_fig,'Pointer','watch');
  
  [this, mask, f] = imroi(D{1});
  close(f);

  % upload new Data sets
  mifit_List_Data_push(mask);

  % apply the mask to all data sets
  for index=1:numel(D)
    this = D{index};
    try
      this = this .* mask;  % can be used with non matching dimensions/axes.
      mifit_List_Data_push(this);
    end
  end
  
  set(mifit_fig,'Pointer','arrow');
