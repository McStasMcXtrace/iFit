function [EXP, fig] = ResLibCal_fig2EXP(fig)
% EXP=ResLibCal_fig2EXP(fig): collect EXP from main ResLibCal window
%
% collect the current data from the figure and build a ResLib EXP structure
%  when not specified, the figure handle is determined from open figures.
%
% Return:
%  EXP: ResLib EXP structure with instrument parameters, or return [] if not opened

% Calls: none

  EXP = [];
  if nargin == 0, fig = ''; end
  if ishandle(fig), 
    if ~strcmp(get(fig, 'Tag'),'ResLibCal'), fig = ''; end
  else fig = ''; 
  end
  if isempty(fig), fig=ResLibCal_fig; end
  if isempty(fig), return; end

  % from ResLib/MakeExp A. Zheludev, 1999-2006 ---------------------------------
  
  %-------------------------   Computation type    -----------------------------
  EXP.method=get(ResLibCal_fig('EXP_method'),'String');
  EXP.method=EXP.method{get(ResLibCal_fig('EXP_method'),'Value')};
  EXP.moncor=0; % Intensity NOT normalized to flux on monitor by default...

  %-------------------------   Monochromator and analyzer    -------------------
  EXP.mono.d       =str2double(get(ResLibCal_fig('EXP_mono_d'),'String'));
  EXP.mono.mosaic  =str2double(get(ResLibCal_fig('EXP_mono_mosaic'),'String'));
  EXP.mono.vmosaic =str2double(get(ResLibCal_fig('EXP_mono_vmosaic'),'String'));
  EXP.ana.d        =str2double(get(ResLibCal_fig('EXP_ana_d'),'String'));
  EXP.ana.mosaic   =str2double(get(ResLibCal_fig('EXP_ana_mosaic'),'String'));
  EXP.ana.vmosaic  =str2double(get(ResLibCal_fig('EXP_ana_vmosaic'),'String'));
  
  if EXP.mono.vmosaic==0, EXP.mono.vmosaic=EXP.mono.mosaic; end
  if EXP.ana.vmosaic==0,  EXP.ana.vmosaic =EXP.ana.mosaic; end
  EXP.mono.tau = 2*pi/EXP.mono.d;
  EXP.ana.tau  = 2*pi/EXP.ana.d;

  %-------------------------   Sample ------------------------------------------
  abc      =str2num(get(ResLibCal_fig('EXP_sample_abc'),'String'));
  if length(abc) == 1, abc = [ abc abc abc ]; end
  if length(abc) == 2, abc = [ abc abc(2)  ]; end
  if length(abc) == 3 && all(isfinite(abc))
    EXP.sample.a = abc(1);
    EXP.sample.b = abc(2);
    EXP.sample.c = abc(3);
    % update value of object
    set(ResLibCal_fig('EXP_sample_abc'),'String',sprintf('%g ',abc));
  end
  abc      =str2num(get(ResLibCal_fig('EXP_sample_alphabetagamma'),'String'));
  if length(abc) == 1, abc = [ abc abc abc ]; end
  if length(abc) == 2, abc = [ abc abc(2)  ]; end
  if length(abc) == 3 && all(isfinite(abc))
    EXP.sample.alpha = abc(1);
    EXP.sample.beta  = abc(2);
    EXP.sample.gamma = abc(3);
    % update value of object
    set(ResLibCal_fig('EXP_sample_alphabetagamma'),'String',sprintf('%g ',abc));
  end
  EXP.sample.mosaic   =str2double(get(ResLibCal_fig('EXP_sample_mosaic'),'String'));
  EXP.sample.vmosaic  =str2double(get(ResLibCal_fig('EXP_sample_vmosaic'),'String'));
  if EXP.sample.vmosaic==0, EXP.sample.vmosaic=EXP.sample.mosaic; end
  
  % current HKLW position (supports vector input for scans)
  EXP.QH=str2num(get(ResLibCal_fig('EXP_QH'),'String'));
  EXP.QK=str2num(get(ResLibCal_fig('EXP_QK'),'String'));
  EXP.QL=str2num(get(ResLibCal_fig('EXP_QL'),'String'));
  EXP.W =str2num(get(ResLibCal_fig('EXP_W' ),'String'));

  %-------------------------   Fixed neutron energy    -------------------------
  EXP.efixed=str2double(get(ResLibCal_fig('EXP_efixed'),'String'));
  EXP.Kfixed=str2double(get(ResLibCal_fig('EXP_Kfixed'),'String'));
  EXP.Lfixed=str2double(get(ResLibCal_fig('EXP_Lfixed'),'String'));
  if get(ResLibCal_fig('EXP_Kf_button'),'Value')
    EXP.infin=-1;      % Fixed final energy: default value
    EXP.fx   =2;
  else
    EXP.infin=1;       % .. negative for fixed incident energy
    EXP.fx   =1;
  end
  % store as well Ki and Kf
  f  = 0.4826;
  EXP.ki=sqrt(EXP.Kfixed^2+(EXP.fx-1)*f*EXP.W); 
  EXP.kf=sqrt(EXP.Kfixed^2-(2-EXP.fx)*f*EXP.W);
  
  %-------------------------   Soller and neutron guide collimation    ---------
  table = get(ResLibCal_fig('EXP_collimators'),'Data');
  EXP.hcol=table(1,:);         % Horizontal collimation: FWHM minutes of arc
  EXP.vcol=table(2,:);         % Vertical collimation: FWHM minutes of arc
  % handle negative values in HCOL/VCOL for guide coating divergence
  % handle case of negative collimations : guide limited
  index = find(EXP.hcol<0);
  if ~isempty(index)
      EXP.hcol(index) = abs(EXP.hcol(index))*0.1*60*(2*pi/EXP.ki)/0.4247/sqrt(3);
  end
  index = find(EXP.vcol<0);
  if ~isempty(index)
      EXP.vcol(index) = abs(EXP.vcol(index))*0.1*60*(2*pi/EXP.ki)/0.4247/sqrt(3);
  end

  %-------------------------   Experimental geometry    ------------------------
  % Value Popup: 1=right -> -1 (clock) : 2=left -> +1
  sm = get(ResLibCal_fig('EXP_mono_dir'),'Value');   if sm==1, sm=-1; else sm=1; end
  ss = get(ResLibCal_fig('EXP_sample_dir'),'Value'); if ss==1, ss=-1; else ss=1; end
  sa = get(ResLibCal_fig('EXP_ana_dir'),'Value');    if sa==1, sa=-1; else sa=1; end
  % scattering direction: 1=opposite to previous element, -1=same as previous
  EXP.mondir = sm;
  if (ss*sm < 0) EXP.dir1 = 1; else EXP.dir1 = -1; end
  if (sm*sa < 0) EXP.dir2 = 1; else EXP.dir2 = -1; end
  EXP.mono.dir  =sm;
  EXP.sample.dir=ss;
  EXP.ana.dir   =sa;
  
  hObject = {'orient1','orient2'};
  for index=1:length(hObject)
    abc      =str2num(get(ResLibCal_fig([ 'EXP_' hObject{index} ]),'String'));
    if length(abc) == 1, 
      if abc, abc = [ abc 0 0 ]; else abc = [ 0 1 0 ]; end
    end
    if length(abc) == 2, abc = [ abc abc(2)  ]; end
    if length(abc) == 3 && all(isfinite(abc))
      EXP.(hObject{index}) = abc;
      % update value of object
      set(ResLibCal_fig([ 'EXP_' hObject{index} ]),'String',num2str(abc));
    end
  end

  %-------------------------   Horizontally focusing analyzer  -----------------
  RAH = str2double(get(ResLibCal_fig('EXP_ana_rh'),'String'));
  if RAH == 0
    EXP.horifoc=-1; %Flat analyzer
  else
    EXP.horifoc=1;  %Horizontally-focused analyzer
  end

  %-------------------------   Spatial parameters     --------------------------
  EXP.beam.width      =str2double(get(ResLibCal_fig('EXP_beam_width'),'String'));
  EXP.beam.height     =str2double(get(ResLibCal_fig('EXP_beam_height'),'String'));
  EXP.detector.width  =str2double(get(ResLibCal_fig('EXP_detector_width'),'String'));
  EXP.detector.height =str2double(get(ResLibCal_fig('EXP_detector_height'),'String'));
  EXP.mono.width      =str2double(get(ResLibCal_fig('EXP_mono_width'),'String'));
  EXP.mono.height     =str2double(get(ResLibCal_fig('EXP_mono_height'),'String'));
  EXP.mono.depth      =str2double(get(ResLibCal_fig('EXP_mono_depth'),'String'));
  EXP.ana.width       =str2double(get(ResLibCal_fig('EXP_ana_width'),'String'));
  EXP.ana.height      =str2double(get(ResLibCal_fig('EXP_ana_height'),'String'));
  EXP.ana.depth       =str2double(get(ResLibCal_fig('EXP_ana_depth'),'String'));
  EXP.monitor.width =2*str2double(get(ResLibCal_fig('EXP_detector_width'),'String'));
  EXP.monitor.heigth=2*str2double(get(ResLibCal_fig('EXP_detector_height'),'String'));

  EXP.sample.width =str2double(get(ResLibCal_fig('EXP_sample_width'),'String'));
  EXP.sample.height=str2double(get(ResLibCal_fig('EXP_sample_height'),'String'));
  EXP.sample.depth =str2double(get(ResLibCal_fig('EXP_sample_depth'),'String'));
  EXP.sample.shape =diag([ EXP.sample.depth EXP.sample.width EXP.sample.height ].^2/12);

  % Spectrometer arms
  EXP.arms=[ table(3,:) table(3,1)*0.7 ];
  
  % limit  collimation/divergences from distances and size of elements (in minutes)
  EXP.hcol(1) = min(EXP.hcol(1), atan2(EXP.mono.width +EXP.beam.width,    EXP.arms(1))*180/2/pi*60);
  EXP.vcol(1) = min(EXP.vcol(1), atan2(EXP.mono.height+EXP.beam.height,   EXP.arms(1))*180/2/pi*60);
  EXP.hcol(2) = min(EXP.hcol(2), atan2(EXP.sample.width +EXP.mono.width,  EXP.arms(2))*180/2/pi*60);
  EXP.vcol(2) = min(EXP.vcol(2), atan2(EXP.sample.height+EXP.mono.height, EXP.arms(2))*180/2/pi*60);
  EXP.hcol(3) = min(EXP.hcol(3), atan2(EXP.ana.width +EXP.sample.width,   EXP.arms(3))*180/2/pi*60);
  EXP.vcol(3) = min(EXP.vcol(3), atan2(EXP.ana.height+EXP.sample.height,  EXP.arms(3))*180/2/pi*60);
  EXP.hcol(4) = min(EXP.hcol(4), atan2(EXP.ana.width +EXP.detector.width, EXP.arms(4))*180/2/pi*60);
  EXP.vcol(4) = min(EXP.vcol(4), atan2(EXP.ana.height+EXP.detector.height,EXP.arms(4))*180/2/pi*60);
  
  % Crystal curvatures [cm-1]
  EXP.mono.rv=str2double(strtok(get(ResLibCal_fig('EXP_mono_rv'),'String')));
  EXP.mono.rh=str2double(strtok(get(ResLibCal_fig('EXP_mono_rh'),'String')));
  EXP.ana.rv =str2double(strtok(get(ResLibCal_fig('EXP_ana_rv'),'String')));
  EXP.ana.rh =str2double(strtok(get(ResLibCal_fig('EXP_ana_rh'),'String')));
  % handle flat mono/ana
  if EXP.mono.rv == 0, EXP.mono.rv=Inf; end
  if EXP.mono.rh == 0, EXP.mono.rh=Inf; end
  if EXP.ana.rv  == 0, EXP.ana.rv=Inf; end
  if EXP.ana.rh  == 0, EXP.ana.rh=Inf; end
  % handle automatic mono/ana curvatures
  if any([ EXP.mono.rv EXP.mono.rh EXP.ana.rv EXP.ana.rh ] < 0)
    rho = rc_focus(EXP);  % in 'private', from ResCal5
    t='';
    if EXP.mono.rv < 0,
      EXP.mono.rv=rho.mh;
      t=[t sprintf('RMV=%g [cm] ', EXP.mono.rv) ];
      set(ResLibCal_fig('EXP_mono_rv'),'String', sprintf('-1 (%g)', EXP.mono.rv));
    end
    if EXP.mono.rh < 0,
      EXP.mono.rh=rho.mv;
      t=[t sprintf('RMH=%g [cm] ', EXP.mono.rh) ];
      set(ResLibCal_fig('EXP_mono_rh'),'String', sprintf('-1 (%g)', EXP.mono.rh));
    end
    if EXP.ana.rv  < 0,
      EXP.ana.rv =rho.ah;
      t=[t sprintf('RAV=%g [cm] ', EXP.ana.rv) ];
      set(ResLibCal_fig('EXP_ana_rv'),'String', sprintf('-1 (%g)', EXP.ana.rh));
    end
    if EXP.ana.rh  < 0,
      EXP.ana.rh =rho.av;
      t=[t sprintf('RAH=%g [cm] ', EXP.ana.rh) ];
      set(ResLibCal_fig('EXP_ana_rh'),'String', sprintf('-1 (%g)', EXP.ana.rv));
    end
    if ~isempty(t)
      disp([ 'Auto focus: ' t]);
    end
  end

% end ResLibCal_fig2EXP
