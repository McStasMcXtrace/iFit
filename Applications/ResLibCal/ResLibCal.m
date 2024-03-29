function out = ResLibCal(varargin)
% ResLibCal: compute and display the neutron TAS resolution function
%
% To start the application and open the main GUI window, use
%   ResLibCal
% To compute directly the resolution function, sending an EXP ResLib-like
%   configuration structure, use:
%   out = ResLibCal(EXP);
% To convolve an iFunc model with a 4D resolution function, use:
%   out = ResLibCal(model);
% To use ResLibCal from the command line, use:
%   out = ReslibCal(command, arguments...);
% where 'command' is one of:
%
%   open    open a configuration file (par, cfg, res, ini, m): ResLibCal('open','file')
%   save    save the configuration in the Preference directory (ini format)
%   saveas  save the configuration into a specified file/format: ResLibCal('saveas','file')
%   export  dump the main ResLibCal window into a file
%   exit    close all active views, and save current configuration
%   reset   re-load the default configuration
%   print   generate an HTML document to be printed
%   create  open the main GUI (start interface), and read last saved configuration
%   compute only compute the matrix (no plotting/printing)
%   update  compute, and then update open views, or send result to the console
%   view2   display the 2D view (resolution projections)
%   view3   display the 3D view (resolution)
%   tas     display the spectrometer geometry
%   default same as create, but does not read last configuration (using reset configuration)
%   quit    same as exit, but does not save the configuration
%   close   close the 2D, 3D and TAS view windows
%   resol   print-out the resolution matrix a la RESCAL
%   bragg   print-out the Bragg widths a la RESCAL
%   list    print-out the RESCAL parameter list
%   config  return the current configuration (ResLib EXP)
%   hkle    return the current HKLE location. Set it back with ResLibCal(hkle{:});
%   silent  use silent computation (no plot/display) for further arguments
%   <PAR>=<VALUE> sets a parameter value, e.g. 'DM=3.355'
%
% To compute the resolution at a given HKLW location, using the current settings
%   resolution = ResLibCal(QH,QK,QL,W)
% where QH,QK,QL,W can be vectors, or empty to use current settings.
% To only compute the resolution function, use:
%   ResLibCal('close'); out=ResLibCal('compute');
% which will close the 2D,3D and TAS views, then compute the resolution
% function, returned in out.resolution
%
% The application contains a main interface with:
% * Menu, Method, Scan and Instrument parameters (main)
% * Resolution function plot (2D)
% * Resolution function plot (3D)
% * Instrument view
%
% when changing any value in the main GUI:
% * Method and Scan parameters, Instrument parameters
% any opened view is updated after a re-computation of the resolution.
%
% The 2D and 3D views can be closed without ending the application.
% When the main window is closed, or Exit is selected all views are closed
%
% Convolution in 4D for TAS ----------------------------------------------------
%
%   The 4D convolution syntax allow to simulate a TAS scan from a parametrised 
% dispersion. In the following example, we simulate a scan through a cubic crystal
% dispersion, and show the ideal S(q,w) as well as the measured, broadened 
% measurement. Then, the simulated scan is inserted in a representation of the 
% S(q,w), to visualise the scan trajectory and simulated signal.
% To use this tool, you need to input a 2D or 4D dispersion model (iFunc).
% The dispersion is either a 2D model S(|q|,w) or a 4D model S(qh,qk,ql,w).
% The axes of the dispersion are in the lattice reciprocal space, in r.l.u.
%
%   s=sqw_cubic_monoatomic; % create a 4D S(q,w) for a cubic pure material
%   t=ResLibCal(s);         % convolute it with a TAS resolution, and open ResLibCal.
%   w=linspace(0.01,20,50); qh=0.3*ones(size(w)); qk=0*qh; ql=qk; % a scan
%   signal1=iData(t, [], qh,qk,ql,w);
%   signal0=iData(s, [], qh,qk,ql,w);
%   figure; plot(squeeze([signal1 signal0*100])); % plot the dispersion and simulated measurement
%   % now plot the 4D dispersion with the scan in it, for fun
%   qh=linspace(0.01,.5,50);qk=qh; ql=qh; w=linspace(0.01,10,51); % a 4D grid
%   f=iData(s,[],qh,qk,ql,w); % evaluate the model on the 4D grid
%   figure; surf(log(f(:,:,1,:)),'half'); hold on;  % plot dispersion, and scan
%   scatter3(log(signal1(:,:,1,:)),'filled');
%
% syntax:
%   ResLibCal('command', arguments...)
%
% input: any combination of:
%   command: string among those listed above, which can be followed by any other
%            allowed parameter.
%   qh,qk,ql,w: 4 vectors or 4D matrices which delimit a region in the reciprocal
%            space where the TAS resolution function should be computed.
%   model:   an iFunc model which is to be convolved with the TAS response.
%   EXP:     a structure holding a ResLibCal, ResCal or ResLib configuration.
% output:
%   a ResLibCal configuration with e.g. 'resolution' field, or a 4D convoluted model.
%
% Example: ResLibCal('Bragg'); isnumeric(r) && numel(r) == 4
%
% Version: $Date$ $Version$ $Author$
% 
% See also iFunc/conv

% Contributions:
% ResCal5: rc_cnmat rc_popma rc_int rc_focus rc_bragg rc_bragghklz
%          rc_re2rc rc_phon
% ResLib:  ResMatS CleanArgs StandardSystem
%          modvec scalar ResMat GetLattice star GetTau
%          element sizes: the sqrt(12) has been removed everywhere. It is taken into account in modified ResMat.

% Private functions: 'out' is a full ResLibCal configuration, 'EXP' is ResLib structure
% out         = ResLibCal_Compute(EXP or out)
% resolution  = ResLibCal_ComputeResMat(EXP or out)
% fig         = ResLibCal_EXP2fig(EXP or out, fig)
% [p, labels] = ResLibCal_EXP2RescalPar(EXP or out)
% [EXP, fig]  = ResLibCal_fig2EXP(fig)
% [res, inst] = ResLibCal_FormatString(out)
% EXP         = ResLibCal_Open(filename, EXP or out) % update EXP/out from file
% EXP         = ResLibCal_RescalPar2EXP(str, EXP)
% RMS         = ResLibCal_RM2RMS(H,K,L,W,EXP,RM)
% EXP         = ResLibCal_SampleRotateS(H,K,L,EXP)
% filename    = ResLibCal_Saveas(filename, EXP)
% ResLibCal_UpdateDTau(handle)
% ResLibCal_UpdateEKLfixed(handle)
% filename = ResLibCal_Save
% out      = ResLibCal_UpdateViews(out)
% out      = ResLibCal_ViewResolution(out, dim)
% out      = ResLibCal_UpdateResolution2(out)
% out      = ResLibCal_UpdateResolution3(out)
% ResLibCal_UpdateTauPopup(handle, EXP)
% ResLibCal_GenerateReport(filename)

out = [];
silent_mode = 0;
ResLibCal_version = [ mfilename ' 1.3 ($Date$)' ];

persistent fig

if nargin == 0
  % nargin == 0

  % open or create the main GUI
  out = feval(mfilename, 'create'); % load last configuration, and Compute
  out = ResLibCal_ViewResolution(out,2);  % open/raise View Res2
  out = ResLibCal_UpdateViews(out); % when they exist
elseif nargin >= 1 && (isa(varargin{1}, 'iFunc') || isa(varargin{1}, 'iData'))
  out = ResLibCal_tas_conv4d(varargin{:});
  return
end
% menu actions:
config = [];
while ~isempty(varargin)
  if isscalar(varargin{1}) && ishandle(varargin{1}) ...
    && (numel(varargin) == 2 && ~isnumeric(varargin{2}))
    % used to update the Table 'collimators'
    varargin = [ {'update_handle'} varargin ];
  end
  if ischar(varargin{1})
    action = varargin{1};
    switch lower(action)
    % menu items ---------------------------------------------------------------
    case {'file_open','open','load'}
      [EXP, fig] = ResLibCal_fig2EXP;
      if isfield(EXP,'EXP'), EXP = EXP.EXP; end
      if length(varargin) <= 1, varargin{2} = ''; end
      if length(varargin) <= 2, varargin{3} = EXP; end
      EXP = ResLibCal_Open(varargin{2:3}, silent_mode);  % (filename, EXP)
      if isempty(config)
        config=ResLibCal_GetConfig(silent_mode);
      end
      out = config;
      if ~isfield(EXP,'EXP')
        % add current EXP to out.EXP;
        out.EXP=mergestruct(out.EXP, EXP);
      else
        % add current EXP to out
        out=mergestruct(out,EXP); 
      end
      out = ResLibCal_Compute(out);
      if ~silent_mode, out = ResLibCal_UpdateViews(out); end
      varargin(2:3) = [];
    case 'file_open_instr'
      p = fileparts(which(mfilename));
      out = ResLibCal('open',[ p filesep 'instruments' ]);
    case 'file_reset'
      filename = fullfile(prefdir, 'ResLibCal.ini');
      if ~exist(filename, 'file')
        source = 'the factory defaults.';
      else
        source = [ 'the file ' filename ' (delete it to return to factory defaults)' ];
      end
      options.Default     = 'Cancel';
      options.Interpreter = 'tex';
      ButtonName = questdlg({ ...
        '{\fontsize{14}{\color{blue}Reload default configuration ?}}', ...
        'This will reset all fields of the ResLibCal window from ', ...
        source, ...
        '{\bf Reset now ?}'}, 'ResLibCal: Reset ?', ...
        'Reset', 'Cancel', 'Factory settings', ...
        options);
      if ~strcmp(ButtonName, 'Cancel')
        if strcmp(ButtonName, 'Reset')
          out = ResLibCal('reset');
        else
          out = ResLibCal('default');
        end
        out = ResLibCal('update'); % compute and update plots plot
      end
    case {'file_save','save'}
      % save configuration so that it is re-opened at next re-start
      filename = ResLibCal_Saveas(fullfile(prefdir, 'ResLibCal.ini'));
    case {'file_saveas','saveas'}
      % save configuration
      ResLibCal_Saveas(out); % (filename, EXP)
      return;
    case {'file_print','print'}
      %fig = ResLibCal_fig;
      %printdlg(fig);
      ResLibCal_GenerateReport;
      return;
    case {'file_export','export'}
      [filename, pathname] = uiputfile( ...
         {'*.pdf',  'Portable Document Format (*.pdf)'; ...
          '*.eps',  'Encapsulated Postscript (*.eps)'; ...
          '*.png',  'Portable Network Graphics image (*.png)'; ...
          '*.jpg',  'JPEG image (*.jpg)'; ...
          '*.tif',  'TIFF image, compressed (*.tif)'; ...
          '*.bmp',  'Windows bitmap (*.bmp)'; ...
          '*.html', 'Full report in Hypertext Markup Language document (*.html)'; ...
          '*.*',  'All Files (*.*)'}, ...
          'ResLibCal: Export configuration window as...');
      if isempty(filename) || all(filename == 0), return; end
      filename = fullfile(pathname, filename);
      [~,~,e] = fileparts(filename);
      if strcmp(e, '.html')
        ResLibCal_GenerateReport(filename);
      else
        fig = ResLibCal_fig;
        saveas(fig, filename);
        disp([ '% Exported ' ResLibCal_version ' main window to file ' filename ]);
      end
    case {'publish','html','report'}
      % create an HTML report
      ResLibCal_GenerateReport;
    case {'file_exit','exit','quit'}
      if isempty(ResLibCal_fig), return; end
      if ~strcmp(action, 'quit')
        % save configuration so that it is re-opened at next re-start
        filename = ResLibCal_Saveas(fullfile(prefdir, 'ResLibCal.ini'));
      end
      % close windows
      hObjects={'ResLibCal',...
                'ResLibCal_View2',...
                'ResLibCal_View3',...
                'ResLibCal_View1'};
      for index=1:length(hObjects)
        fig=findall(0, 'Tag',hObjects{index});
        if ishandle(fig), delete(fig); end
      end
    case {'view_resolution2','view2'}
      if length(varargin) > 1, v=varargin{2}; varargin(2) = []; else v=[]; end
      out = ResLibCal_ViewResolution(v,2);  % open/raise View Res2
      out = ResLibCal_UpdateViews(out, 'force');
    case {'view_resolution3','view3'}
      if length(varargin) > 1, v=varargin{2}; varargin(2) = []; else v=[]; end
      out = ResLibCal_ViewResolution(v,3);  % open/raise View Res3
      out = ResLibCal_UpdateViews(out, 'force');
    case {'view_tas','geometry','tas'}
      if length(varargin) > 1, v=varargin{2}; varargin(2) = []; else v=[]; end
      out = ResLibCal_ViewResolution(v,1);  % open/raise View TAS
      out = ResLibCal_UpdateViews(out, 'force');
    case {'help_content','help'}
      link = fullfile(fileparts(which(mfilename)), 'doc', [ mfilename '.html' ]);
      disp([ mfilename ': opening help from ' link ])
      web(link);
      out = link;
      return;
    case 'version'
      message = [ ResLibCal_version ' compute and display the triple-axis ' ...
        'resolution function obtained from e.g. Cooper-Nathans and Popovici ' ...
        'analytical approximations. Part of <ifit.mccode.org>. ' ...
        '(c) E. Farhi, ILL, Soleil. EUPL. ' ...
        'Contributions from: A. Zheludev, A. Tennant, D. Mc Morrow, J. Ollivier, B. Hennion, Hargreave,Hullah, N. Moshtagh' ];

      out = message;
      return;
    case {'help_about'}
      % get the ILL logo from object
      fig = ResLibCal_fig;
      
      if isempty(fig), return; end
      cdata = get(ResLibCal_fig('ILL_logo'), 'CData');
      message = {...
        [ '{\fontsize{14}{\color{blue}' ResLibCal_version '}}' ], ...
        'ResLibCal is a graphical user interface to compute and display' , ...
        'the {\bf triple-axis resolution function} obtained from e.g. Cooper-Nathans and Popovici analytical approximations. The GUI allows to select among a set of computation kernels.' , ...
        'This application was written by E. Farhi {\copyright}ILL, Soleil <emmanuel.farhi@synchrotron-soleil.fr> (EUPL license) using' , ...
        '\bullet ResLib 3.4 (A. Zheludev)' , ...
        '\bullet ResCal5 (A. Tennant and D. Mc Morrow)' , ...
        '\bullet Res3ax (J. Ollivier)' , ...
        '\bullet Rescal/AFILL (Hargreave, Hullah) and vTAS/PKFIT view (A. Bouvet/A. Filhol)', ...
        '\bullet McStas (K. Nielsen, K. Lefmann, E. Farhi, P. Willendrup, et al)', ...
        'Additional contributions from: B. Hennion, N. Moshtagh' };
      CreateMode.WindowStyle = 'modal';
      CreateMode.Interpreter='tex';
      msgbox(message, ...
        'About: ResLibCal', ...
        'custom',cdata,jet(64),CreateMode);
    case 'view_update'
      if numel(varargin) > 1 && isstruct(varargin{2})
        out = ResLibCal_Compute(varargin{2}); % arg can be an EXP
        varargin(2) = [];
      else
        out = ResLibCal_Compute;
      end
      ResLibCal_ViewResolution(out,2); % if not opened, open at least the 2D view
      ResLibCal_UpdateViews(out, 'force');
    case {'view_autoupdate','autoupdate'}
      status = '';
      if numel(varargin) > 1 && ischar(varargin{2})
        status = varargin{2};
        if strcmp(status,'on') || strcmp(status,'off')
          varargin(2) = [];
        else
          status = ''; % use current setting and toggle
        end
      end
      if isempty(status)  % toggle
        status = get(ResLibCal_fig('View_AutoUpdate'), 'Checked');
        if strcmp(status,'off'), status = 'on'; end
      end
      if ~strcmp(status, 'on'), status = 'off'; end % make sure we get on or off
      set(ResLibCal_fig('View_AutoUpdate'), 'Checked', status);
    case {'view_resolutionrlu','rlu'}
      set(ResLibCal_fig('View_ResolutionRLU'), 'Checked', 'on');
      set(ResLibCal_fig('View_ResolutionSPEC'),'Checked', 'off');
      set(ResLibCal_fig('View_ResolutionABC'), 'Checked', 'off');
      set(ResLibCal_fig('View_ResolutionLattice'), 'Checked', 'off');
      ResLibCal_UpdateViews([],'force');
    case {'view_resolutionabc','abc'}
      set(ResLibCal_fig('View_ResolutionRLU'), 'Checked', 'off');
      set(ResLibCal_fig('View_ResolutionSPEC'),'Checked', 'off');
      set(ResLibCal_fig('View_ResolutionABC'), 'Checked', 'on');
      set(ResLibCal_fig('View_ResolutionLattice'), 'Checked', 'off');
      ResLibCal_UpdateViews([],'force');
    case {'view_resolutionspec','spec'}
      set(ResLibCal_fig('View_ResolutionRLU'), 'Checked', 'off');
      set(ResLibCal_fig('View_ResolutionSPEC'),'Checked', 'on');
      set(ResLibCal_fig('View_ResolutionABC'), 'Checked', 'off');
      set(ResLibCal_fig('View_ResolutionLattice'), 'Checked', 'off');
      ResLibCal_UpdateViews([],'force');
    case {'view_resolutionlattice','lattice'}
      set(ResLibCal_fig('View_ResolutionRLU'), 'Checked', 'off');
      set(ResLibCal_fig('View_ResolutionSPEC'),'Checked', 'off');
      set(ResLibCal_fig('View_ResolutionABC'), 'Checked', 'off');
      set(ResLibCal_fig('View_ResolutionLattice'), 'Checked', 'on');
      ResLibCal_UpdateViews([],'force');
    case {'view_resolutionhv','hv'}
      % select the 2nd axis. First is always a* or A (1)
      % second can be selected here as index 2 or 3 in 2D 3D view.
      % then the resolutionxyz=Q below should select the other one
      status = get(ResLibCal_fig('View_ResolutionHV'), 'Checked');
      if strcmp(status,'on'), 
        status = 'off'; 
        set(ResLibCal_fig('View_ResolutionHV'), 'Label','Resolution: horizontal in [Qx,Qy]');
      else 
        status = 'on'; 
        set(ResLibCal_fig('View_ResolutionHV'), 'Label','Resolution: vertical in [Qx,Qz]');
      end
      set(ResLibCal_fig('View_ResolutionHV'), 'Checked', status);
      ResLibCal_UpdateViews([],'force');
    case {'view_resolutionxyz','zw'}
      status = get(ResLibCal_fig('View_ResolutionXYZ'), 'Checked');
      if strcmp(status,'on'), 
        status = 'off'; 
        set(ResLibCal_fig('View_ResolutionXYZ'), 'Label','Resolution: vertical in [E]');
      else 
        status = 'on'; 
        set(ResLibCal_fig('View_ResolutionXYZ'), 'Label','Resolution: vertical in [Q]');
      end
      set(ResLibCal_fig('View_ResolutionXYZ'), 'Checked', status);
      ResLibCal_UpdateViews([],'force');
    case 'view_resolution_cloud'
      status = get(ResLibCal_fig('View_Resolution_Cloud'), 'Checked');
      if strcmp(status,'on'), status = 'off'; else status = 'on'; end
      set(ResLibCal_fig('View_Resolution_Cloud'), 'Checked', status);
      ResLibCal_UpdateViews([],'force');
    case {'view_nmc','nmc','monte-carlo'}
      if isfield(out,'EXP'), EXP=out.EXP; else EXP=out; end
      if     isfield(EXP, 'NMC'),         NMC=EXP.NMC;
      else
        NMC=get(ResLibCal_fig('View_NMC'), 'UserData');
      end
      if isempty(NMC) || ~isnumeric(NMC), NMC  = 200; end
      NMC = inputdlg('Enter the number of Monte-Carlo iterations (cloud)', ...
        'ResLibCal: Monte-Carlo iterations ?',1,{ num2str(NMC) });
      if ~isempty(NMC)
        NMC=str2double(NMC{1});
        if NMC > 0
          if isfield(out,'EXP'), out.EXP.NMC=NMC; else EXP.NMC=NMC; end
          set(ResLibCal_fig('View_NMC'), 'UserData', NMC);
        end
      end
    case {'view_close','close'}
      delete(findobj(0, 'Tag', 'ResLibCal_View1'));
      delete(findobj(0, 'Tag', 'ResLibCal_View2'));
      delete(findobj(0, 'Tag', 'ResLibCal_View3'));
      
    % RESCAL actions -----------------------------------------------------------
    case 'list'
      out = ResLibCal_Compute(out);
      disp('RESCAL parameters')
      disp(out.ResCal);
    case 'bragg'
      out = ResLibCal_Compute(out);
      resolution = out.resolution;
      if isstruct(resolution), resolution = { resolution }; end
      for index=1:numel(resolution)
        H   = resolution{index}.HKLE(1); K=resolution{index}.HKLE(2); 
        L   = resolution{index}.HKLE(3); W=resolution{index}.HKLE(4);
        fprintf(1,'QH=%5.3g QK=%5.3g QL=%5.3g [rlu] E=%5.3g [meV]\n', H,K,L,W);
        disp('  BRAGG Widths, Radial,tangential, Vertical (HWHM) [ANG-1]');
        FrameStr = resolution{index}.spec.frameStr;
        disp(['  Rad:along A=' FrameStr{1} '; Tang:along ' FrameStr{2} '; Vert:along ' FrameStr{3} ])
        fprintf(1, '  DQR=%g DQT=%g DQV=%g\n', resolution{index}.spec.Bragg(1:3)/2);
        disp('  Energy Widths (HWHM) [meV]');
        fprintf(1, '  DVN=%g DEE=%g\n', resolution{index}.spec.Bragg([ 5 4 ])/2);
        disp('----------------------------------------------------------');
        out = resolution{index}.spec.Bragg;
      end
    case 'resol'
      out = ResLibCal_Compute(out);
      resolution = out.resolution;
      if numel(out.resolution) == 1, resolution={ out.resolution }; end
      for index=1:numel(resolution)
        H   = resolution{index}.HKLE(1); K=resolution{index}.HKLE(2); 
        L   = resolution{index}.HKLE(3); W=resolution{index}.HKLE(4);
        fprintf(1,'QH=%5.3g QK=%5.3g QL=%5.3g [rlu] E=%5.3g [meV]\n', H,K,L,W);
        disp('----------------------------------------------------------');
        % display the resolution matrix in all available frames
        for frames={'rlu','spec','ABC'}  % others: 'cart','rlu_ABC','ABC'
          frame = resolution{index}.(frames{1});
          disp(' ');
          disp([ 'Resolution Matrix [' frames{1} '] ' frame.README ]);
          disp([ 'X:' frame.frameStr{1} '; Y:' frame.frameStr{2} '; Z:' frame.frameStr{3} ]);
          disp('    X        Y        Z        W')
          disp(num2str(frame.RM,'%.1f '));
        end
        disp('----------------------------------------------------------');
      end
    % other actions (not menu items) -------------------------------------------
    case 'silent'
      silent_mode = 1;
    case 'default'  % factory default
      fig = ResLibCal_fig;
      if ~isempty(fig) && ishandle(fig)
        delete(fig);
      end
      f=openfig('ResLibCal');
      out = ResLibCal_Compute;
      % close figure again if it was not there (pure batch mode)
      if isempty(fig) && numel(varargin) == 0, delete(f); end
    case 'reset'    % restore settings from ini file (when exists) or default
      filename = fullfile(prefdir, 'ResLibCal.ini');
      if exist(filename, 'file')
        out = ResLibCal_Open(filename, [], silent_mode); % open the 'ResLibCal.ini' file (last saved configuration)
        out = ResLibCal_Compute(out);
      else
        out = ResLibCal('default');
      end
      out = ResLibCal_UpdateViews(out);
    case 'create' % open GUI if not there yet with config file
      fig = ResLibCal_fig;
      if isempty(fig) || ~ishandle(fig)
        disp([ 'Welcome to ' ResLibCal_version ]);
        fig = openfig('ResLibCal'); % open the main ResLibCal figure.
        set(fig, 'NextPlot','new'); % protect from plotting on top
        if strcmp(action, 'create') % default: ignore config file
          filename = fullfile(prefdir, 'ResLibCal.ini');
          out = ResLibCal_Open(filename); % open the 'ResLibCal.ini' file (last saved configuration)
        end
      elseif length(fig) > 1
        delete(fig(2:end)); % remove duplicated windows
      end
      out = ResLibCal_Compute;
    case 'update' % (this is called when changing the computational method in the GUI)
      % update all opened views with new computation (widget update)
      
      out = ResLibCal_Compute;
      
      if ~silent_mode
        fig = ResLibCal_fig;
        if ~isempty(fig)
          out = ResLibCal_UpdateViews(out);
        elseif nargout == 0
          out = ResLibCal_UpdateViews(out, 'stdout'); % display result to stdout
        end
      end
    case {'compute','resolution'}
      % only compute. No output except in varargout
      % fig = ResLibCal_fig;
      % if no interface exists, load the last saved configuration before computing
      out = ResLibCal_Compute(out);
    case 'update_d_tau'
      % update d-spacing from a popup item
      ResLibCal_UpdateDTau(varargin{2});      % arg is popup handle
      varargin(2) = [];
    case 'update_d_tau_popup'
      % update the mono-ana popup from the DM/DA when value is close
      ResLibCal_UpdateTauPopup;
    case 'update_ekl'
      % update E, K, lambda
      ResLibCal_UpdateEKLfixed(varargin{2});  % arg is edit box handle
      varargin(2) = [];
    case 'update_handle'
      % handle uitable/popup/menu and other updates from widget CallBack
      if length(varargin) < 2, return; end % requires handle as arg
      h     = varargin{2};
      tag   = get(h, 'Tag');
      % handle case of uitable collimators/distances
      if strcmp(tag,'EXP_collimators') && strcmp(get(h,'Type'), 'uitable')
        event = varargin{3};
        data  = get(h,'Data'); % NewData
        if ~isempty(event.Error) || isnan(event.NewData)
          % revert invalid to previous data
          data(event.Indices(1),event.Indices(2)) = event.PreviousData;
          set(h,'Data',data);
        end
        varargin(3)=[];
      elseif any(strcmp(tag,{'EXP_mono_tau_popup','EXP_ana_tau_popup'}))
        ResLibCal_UpdateDTau(h);
      elseif any(strcmp(tag,{'EXP_efixed','EXP_Kfixed','EXP_Lfixed'}))
        ResLibCal_UpdateEKLfixed(h);
      elseif any(strcmp(tag,{'EXP_mono_d','EXP_ana_d'}))
        ResLibCal_UpdateTauPopup;
      end
      varargin(2)=[];
      % update computation and plots
      feval(mfilename, 'update');
    case {'config','EXP'}
      if isempty(config)
        config=ResLibCal_GetConfig(silent_mode);
      end
      out = config;
    case 'hkle'
      if ~isempty(ResLibCal_fig)
        out = { str2num(get(ResLibCal_fig('EXP_QH'),'String'))
          str2num(get(ResLibCal_fig('EXP_QK'),'String'))
          str2num(get(ResLibCal_fig('EXP_QL'),'String'))
          str2num(get(ResLibCal_fig('EXP_W'),'String')) };
      else
        out = {};
      end
      return;
    case 'compile'
      try
        out = ResLibCal_RM2clouds_mcstas('compile');
      catch ME
        out = 'FAILED';
        disp(getReport(ME));
      end
      disp([ mfilename ': using McStas/templateTAS executable: ' out ]);
      return
    otherwise
      % open file name or list of parameters given as 'VAR=VAL; ...'  
      if numel(varargin) > 1 && isstruct(varargin{2})
        config = varargin{2};
        out = config;
        varargin(2)=[];
      end
      if isempty(out), 
        if isempty(config)
          config=ResLibCal_GetConfig(silent_mode);
        end
        out = config;
      end
      if ~isempty(action)
        out = ResLibCal_Open(action, out, silent_mode); % update 'out/EXP' from file
        ResLibCal_EXP2fig(out);                        % put it into the main GUI
      elseif isempty(out)
        out = ResLibCal_GetConfig(silent_mode);
      end
      if numel(varargin) == 0
        out = ResLibCal_Compute(out); % compute the resolution
        if ~silent_mode, 
          ResLibCal_UpdateViews(out); % update views when they exist
        end
      end
      
    end % switch (action as a char command)
    
    varargin(1) = [];
    % end if varargin is char
  elseif isstruct(varargin{1})
    % read an out or EXP structure
    config = varargin{1};
    out = config;
    if ~isfield(out,'EXP')
      EXP=out; out=[];
      out.EXP = EXP; 
    end
    varargin(1) = [];
  elseif isnumeric(varargin{1}) && numel(varargin{1}) == 4 && numel(varargin) == 1
    % ResLibCal([qh qk ql w])
    hkle = varargin{1};
    varargin = { hkle(1) hkle(2) hkle(3) hkle(4) };
    
  elseif numel(varargin) >= 4 && isnumeric(varargin{1}) && isnumeric(varargin{2}) ...
    && isnumeric(varargin{3}) && isnumeric(varargin{4})
    % ResLibCal(qh, qk, ql, w)
    % read HKLE coordinates and compute resolution there
    % get current config
    if isempty(config)
      config = ResLibCal_GetConfig;
    end
    out = config;
    if isfield(out,'EXP') EXP = out.EXP; else 
      if isfield(out,'method') && isfield(out, 'Kfixed'), EXP=out;
      else EXP=[]; end
    end
    if isempty(EXP) || ~isstruct(EXP), return; end

    if ~isempty(varargin) && isnumeric(varargin{1}) 
      if ~isempty(varargin{1}),
        EXP.QH = varargin{1};
        if ~silent_mode,set(ResLibCal_fig('EXP_QH'),'String', mat2str(EXP.QH));end
      end
      varargin(1)=[];
    end
    if ~isempty(varargin) && isnumeric(varargin{1}) 
      if ~isempty(varargin{1}),
        EXP.QK = varargin{1};;
        if ~silent_mode,set(ResLibCal_fig('EXP_QK'),'String', mat2str(EXP.QK));end
      end
      varargin(1)=[];
    end
    if ~isempty(varargin) && isnumeric(varargin{1}) 
      if ~isempty(varargin{1}),
        EXP.QL = varargin{1};
        if ~silent_mode,set(ResLibCal_fig('EXP_QL'),'String', mat2str(EXP.QL));end
      end
      varargin(1)=[];
    end
    if ~isempty(varargin) && isnumeric(varargin{1}) 
      if ~isempty(varargin{1}),
        EXP.W = varargin{1};
        if ~silent_mode,set(ResLibCal_fig('EXP_W'),'String', mat2str(EXP.W));end
      end
      varargin(1)=[];
    end
    out.EXP=EXP;
    if isempty(varargin)
      out = ResLibCal_Compute(EXP);
      if ~silent_mode, 
        if ~isempty(fig), ResLibCal_UpdateViews(out); % when they exist
        elseif nargout==0, ResLibCal_UpdateViews(out,'stdout'); end
      end
    end
  elseif numel(varargin) >= 1 && isempty(varargin{1})
    if nargout, config = ResLibCal_GetConfig(silent_mode); out=config; end
    varargin(1)=[];
  else
    disp([ mfilename ': unknown parameter of class ' class(varargin{1}) ' . Skipping.' ]);
    disp(varargin{1});
    varargin(1)=[];
  end % if type(varargin)
  
end % end while nargin > 0
% end ResLibCal main

% ==============================================================================
function ResLibCal_UpdateTauPopup
% update the popup menu from the editable mono/ana value when d is close
%
  [EXP,fig] = ResLibCal_fig2EXP;
  if isempty(fig), return; end
  popup = ResLibCal_fig('EXP_mono_tau_popup');
  label = GetTau(EXP.mono.tau, 'getlabel');
  index = find(strncmpi(get(popup,'String'), label, length(label)));
  if ~isempty(index) && ~isempty(label), set(popup, 'value', index(1)); end

  popup = ResLibCal_fig('EXP_ana_tau_popup');
  label = GetTau(EXP.ana.tau, 'getlabel');
  index = find(strncmpi(get(popup,'String'), label, length(label)));
  if ~isempty(index) && ~isempty(label), set(popup, 'value', index(1)); end

% ==============================================================================

  
