function mifit_File_Exit(varargin)
% File/Exit: Quit

  mifit_disp([ '[Exit] Exiting miFit. Bye bye.' ]);
  if ishandle(mifit_fig('mifit_View_Parameters'))
    delete(mifit_fig('mifit_View_Parameters'))
  end

  delete(mifit_fig); %  and let Matlab clear all associated memory
  
