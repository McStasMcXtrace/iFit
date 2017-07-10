function ResLibCal_MethodEnableControls(out)
% ResLibCal_MethodEnableControls: activates/disactivates controls from CN to Popovici
%
  fig = ResLibCal_fig;
  if isempty(fig), return; end
  
  Popovici = {'EXP_beam_width', 'EXP_beam_height', ...
  'EXP_detector_width', 'EXP_detector_height', ...
  'EXP_mono_width', 'EXP_mono_height', 'EXP_mono_depth', ...
  'EXP_ana_width', 'EXP_ana_height', 'EXP_ana_depth', ...
  'EXP_sample_width',  'EXP_sample_depth', 'EXP_sample_height', ...
  'EXP_mono_rv', 'EXP_mono_rh', 'EXP_ana_rv', 'EXP_ana_rh'};
  for tag=Popovici
    hObject = ResLibCal_fig(tag{1});
    if ~isempty(hObject)
      if ~isempty(strfind(lower(out.EXP.method), 'popovici')) || ...
         ~isempty(strfind(lower(out.EXP.method), 'mcstas'))
       % this is Popovici or Mcstas method
        set(hObject, 'Enable','on');
      else
        set(hObject, 'Enable','off');
      end
    end
  end
  % special case for Cooper-Nathans legacy without vertical mosaic components
  if ~isempty(strfind(lower(out.EXP.method), 'cooper')) && ...
    (~isempty(strfind(lower(out.EXP.method), 'afill')) || ~isempty(strfind(lower(out.EXP.method), 'rescal5')))
    set(ResLibCal_fig('EXP_mono_vmosaic'), 'Enable','off');
    set(ResLibCal_fig('EXP_ana_vmosaic'), 'Enable','off');
    set(ResLibCal_fig('EXP_sample_vmosaic'), 'Enable','off');
  else
    set(ResLibCal_fig('EXP_mono_vmosaic'), 'Enable','on');
    set(ResLibCal_fig('EXP_ana_vmosaic'), 'Enable','on');
    set(ResLibCal_fig('EXP_sample_vmosaic'), 'Enable','on');
  end
