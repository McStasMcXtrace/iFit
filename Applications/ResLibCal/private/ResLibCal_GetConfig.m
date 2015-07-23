function [out, fig] = ResLibCal_GetConfig(out)
  % get the current configuration or from default file
  % return an out structure with EXP.
  fig = ResLibCal_fig; out=[];
  if ~isempty(fig)
    EXP = ResLibCal_fig2EXP(fig);
    out.Title  = 'ResLibCal';
    out.EXP    = EXP;
    out.handle = fig;
  else
    filename = fullfile(prefdir, 'ResLibCal.ini');
    out = ResLibCal_Open(filename); % open the 'ResLibCal.ini' file (last saved configuration)
    if ~isfield(out,'EXP'),
      EXP = out;
      out.Title  = 'ResLibCal';
      out.EXP    = EXP;
      out.handle = fig;
    end
  end
  
