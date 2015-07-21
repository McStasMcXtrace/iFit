function out = ResLibCal(varargin)
% ResLibCal: compute and display the neutron TAS resolution function
%
% To start the application and open the main GUI window, use
%   ResLibCal
% To compute directly the resolution function, sending an EXP ResLib-like
%   configuration structure, use:
%   out = ResLibCal(EXP);
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
%   print   print the main ResLibCal window
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
% any opened view is updated after a recomputation of the resolution
%
% The 2D and 3D views can be closed without ending the application.
% When the main window is closed, or Exit is selected all views are closed
%
% Version: $Date$

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
%
% Private (inline) functions:
% filename = ResLibCal_Save
% out      = ResLibCal_UpdateViews(out)
% out      = ResLibCal_ViewResolution(out, dim)
% out      = ResLibCal_UpdateResolution2(out)
% out      = ResLibCal_UpdateResolution3(out)
% ResLibCal_UpdateTauPopup(handle, EXP)

out = [];
ResLibCal_version = [ mfilename ' ($Date$)' ];

persistent fig

if nargin == 0
  % nargin == 0

  % open or create the main GUI
  feval(mfilename, 'create'); % load last configuration
  out = ResLibCal_Compute;
  out = ResLibCal_ViewResolution(out,2);  % open/raise View Res2
  out = ResLibCal_UpdateViews(out); % when they exist
end
% menu actions:
while ~isempty(varargin)
  if ishandle(varargin{1})
    varargin = [ {'update_handle'} varargin ];
  end
  if ischar(varargin{1})
    % check if the application window exists, else open it
%    fig=findall(0, 'Tag','ResLibCal');
%    if length(fig) > 1
%      delete(fig(2:end)); % remove duplicated windows
%      fig=fig(1);
%    end
%    if isempty(fig) || ~ishandle(fig)
%      fig = openfig('ResLibCal'); % open the main ResLibCal figure.
%      feval(mfilename, 'create'); % load last configuration
%    end

    action = varargin{1};
    switch lower(action)
    % menu items ---------------------------------------------------------------
    case {'file_open','open','load'}
      [EXP, fig] = ResLibCal_fig2EXP;
      if isfield(EXP,'EXP'), EXP = EXP.EXP; end
      if length(varargin) <= 1, varargin{2} = ''; end
      if length(varargin) <= 2, varargin{3} = EXP; end
      out = ResLibCal_Open(varargin{2:3});  % (filename, EXP)
      out = ResLibCal_Compute(out);
      out = ResLibCal_UpdateViews(out);
      varargin(2:3) = [];
    case 'file_open_instr'
      p = fileparts(which(mfilename));
      ResLibCal('open',[ p filesep 'instruments' ]);
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
      ResLibCal_Save; % (filename=prefdir)
    case {'file_saveas','saveas'}
      % save configuration
      ResLibCal_Saveas(out); % (filename, EXP)
    case {'file_print','print'}
      fig = ResLibCal_fig;
      printdlg(fig);
    case {'file_export','export'}
      [filename, pathname] = uiputfile( ...
         {'*.pdf',  'Portable Document Format (*.pdf)'; ...
          '*.eps',  'Encapsulated Postscript (*.eps)'; ...
          '*.png',  'Portable Network Graphics image (*.png)'; ...
          '*.jpg',  'JPEG image (*.jpg)'; ...
          '*.tif',  'TIFF image, compressed (*.tif)'; ...
          '*.bmp',  'Windows bitmap (*.bmp)'; ...
          '*.*',  'All Files (*.*)'}, ...
          'Export configuration window as...');
      if isempty(filename) || all(filename == 0), return; end
      filename = fullfile(pathname, filename);
      fig = ResLibCal_fig;
      saveas(fig, filename);
      disp([ '% Exported ' ResLibCal_version ' window to file ' filename ]);
    case {'file_exit','exit','quit'}
      if ~strcmp(action, 'quit')
        % save configuration so that it is re-opened at next re-start
        ResLibCal_Save;
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
      out = ResLibCal_UpdateViews(out);
    case {'view_resolution3','view3'}
      if length(varargin) > 1, v=varargin{2}; varargin(2) = []; else v=[]; end
      out = ResLibCal_ViewResolution(v,3);  % open/raise View Res2
      out = ResLibCal_UpdateViews(out);
    case {'view_tas','geometry','tas'}
      if length(varargin) > 1, v=varargin{2}; varargin(2) = []; else v=[]; end
      out = ResLibCal_ViewResolution(v,1);  % open/raise View TAS
      out = ResLibCal_UpdateViews(out);
    case {'help_content','help'}
      link = fullfile(fileparts(which(mfilename)), 'doc', [ mfilename '.html' ]);
      disp([ mfilename ': opening help from ' link ])
      web(link);
      out = link;
    case 'version'
      message = [ ResLibCal_version ' compute and display the triple-axis ' ...
        'resolution function obtained from e.g. Cooper-Nathans and Popovici ' ...
        'analytical approximations. Part of <ifit.mccode.org>.' ...
        'E. Farhi, ILL/Computing for Science.' ];

      out = message;
    case {'help_about'}
      % get the ILL logo from object
      fig = ResLibCal_fig;
      if isempty(fig), return; end
      cdata = get(ResLibCal_fig('ILL_logo'), 'CData');
      message = {...
        [ '{\fontsize{14}{\color{blue}' ResLibCal_version '}}' ], ...
        'ResLibCal is a graphical user interface to compute and display' , ...
        'the {\bf triple-axis resolution function} obtained from e.g. Cooper-Nathans and Popovici analytical approximations. The GUI allows to select among a set of computation kernels.' , ...
        'This application was written by E. Farhi {\copyright}ILL/DS/CS <farhi@ill.eu> (EUPL license) using' , ...
        '\bullet ResLib 3.4 (A. Zheludev)' , ...
        '\bullet ResCal5 (A. Tennant and D. Mc Morrow)' , ...
        '\bullet Res3ax (J. Ollivier)' , ...
        '\bullet Rescal/AFILL (Hargreave, Hullah) and vTAS/PKFIT view (A. Bouvet/A. Filhol)' };
      CreateMode.WindowStyle = 'modal';
      CreateMode.Interpreter='tex';
      msgbox(message, ...
        'About: ResLibCal', ...
        'custom',cdata,jet,CreateMode);
    case 'view_update'
      if numel(varargin) > 1 && isstruct(varargin{2})
        out = ResLibCal_Compute(varargin{2}); % arg can be an EXP
        varargin(2) = [];
      else
        out = ResLibCal_Compute;
      end
      ResLibCal_ViewResolution(out,2); % if not opened, open at least the 2D view
      ResLibCal_UpdateViews(out);
    case 'view_autoupdate'
      status = get(gcbo, 'Checked');
      if strcmp(status,'on'), status = 'off'; else status = 'on'; end
      set(gcbo, 'Checked', status);
    case 'view_resolutionrlu'
      status = get(gcbo, 'Checked');
      if strcmp(status,'on'), status = 'off'; else status = 'on'; end
      set(gcbo, 'Checked', status);
      ResLibCal_UpdateViews;
    case 'view_resolutionxyz'
      status = get(gcbo, 'Checked');
      if strcmp(status,'on'), status = 'off'; else status = 'on'; end
      if strcmp(status,'on')
        set(gcbo, 'Label','Resolution in [Qx,Qy,Qz]');
      else
        set(gcbo, 'Label','Resolution in [Qx,Qy,E]');
      end
      set(gcbo, 'Checked', status);
      ResLibCal_UpdateViews;
    case {'view_close','close'}
      delete(findobj(0, 'Tag', 'ResLibCal_View1'));
      delete(findobj(0, 'Tag', 'ResLibCal_View2'));
      delete(findobj(0, 'Tag', 'ResLibCal_View3'));
      
    % RESCAL actions
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
        FrameStr = strrep(resolution{index}.xyz.FrameStr,'\surd', 'sqrt');
        disp(['  Rad:along A=' FrameStr{1} '; Tang:along ' FrameStr{2} '; Vert:along ' FrameStr{3} ])
        fprintf(1, '  DQR=%g DQT=%g DQV=%g\n', resolution{index}.xyz.Bragg(1:3)/2);
        disp('  Energy Widths (HWHM) [meV]');
        fprintf(1, '  DVN=%g DEE=%g\n', resolution{index}.xyz.Bragg([ 5 4 ])/2);
        disp('----------------------------------------------------------');
      end
    case 'resol'
      out = ResLibCal_Compute(out);
      resolution = out.resolution;
      if numel(out.resolution) == 1, resolution={ out.resolution }; end
      for index=1:numel(resolution)
        H   = resolution{index}.HKLE(1); K=resolution{index}.HKLE(2); 
        L   = resolution{index}.HKLE(3); W=resolution{index}.HKLE(4);
        fprintf(1,'QH=%5.3g QK=%5.3g QL=%5.3g [rlu] E=%5.3g [meV]\n', H,K,L,W);
        disp('  Resolution Matrix, X-AXIS Along Q [ANGS-1] & [meV]')
        FrameStr = strrep(resolution{index}.xyz.FrameStr,'\surd', 'sqrt');
        disp(['  X:along A=' FrameStr{1} '; Y:along ' FrameStr{2} '; Z:along ' FrameStr{3} ])
        disp(' ')
        disp('    X        Y        Z        W')
        disp(num2str(resolution{index}.xyz.RM,'%.1f '));
        disp(' ');
        disp('  Resolution Matrix, Axes WRT Recip. Lattice [R.l.u.] & [meV]')
        FrameStr = strrep(resolution{index}.abc.FrameStr,'\surd', 'sqrt');
        disp(['  X:along A=' FrameStr{1} '; Y:along ' FrameStr{2} '; Z:along ' FrameStr{3} ])
        disp(' ')
        disp('    X        Y        Z        W')
        disp(num2str(resolution{index}.abc.RM,'%.1f '));
        disp('----------------------------------------------------------');
      end
    % other actions (not menu items) -------------------------------------------
    case 'default'  % factory default
      fig = ResLibCal_fig;
      if ~isempty(fig) && ishandle(fig)
        delete(fig);
      end
      f=openfig('ResLibCal');
      out = ResLibCal_Compute;
      % close figure again if it was not there (pure batch mode)
      if isempty(fig), delete(f); end
    case 'reset'    % restore settings from ini file (when exists) or default
      filename = fullfile(prefdir, 'ResLibCal.ini');
      if exist(filename, 'file')
        out = ResLibCal_Open(filename); % open the 'ResLibCal.ini' file (last saved configuration)
        out = ResLibCal_Compute(out);
      else
        out = ResLibCal('default');
      end
      out = ResLibCal_UpdateViews(out);
    case 'create' % open GUI if not there yet with config file
      fig = ResLibCal_fig;
      if isempty(fig) || ~ishandle(fig)
        disp([ 'Welcome to ' ResLibCal_version ]);
        openfig('ResLibCal'); % open the main ResLibCal figure.
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
      fig = ResLibCal_fig;
      
      out = ResLibCal_Compute;
      
      if ~isempty(fig) && strcmp(get(ResLibCal_fig('View_AutoUpdate'), 'Checked'), 'on')
        out = ResLibCal_UpdateViews(out);
      elseif isempty(fig)
        out = ResLibCal_UpdateViews(out); % display result to stdout
      end
    case {'compute','resolution'}
      % only compute. No output except in varargout
      fig = ResLibCal_fig;
      % if no interface exists, load the last saved configuration before computing
      if isempty(fig)
        if isempty(out)
          filename = fullfile(prefdir, 'ResLibCal.ini');
          out = ResLibCal_Open(filename); % open the 'ResLibCal.ini' file (last saved configuration)
        end
        out = ResLibCal_Compute(out);
      else
        out = ResLibCal_Compute;  % get config from main GUI

        ResLibCal_UpdateViews(out); % when they exist
      end
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
        else
          out = feval(mfilename, 'update');
        end
        varargin(3)=[];
      elseif strcmp(get(h,'Type'), 'uitable')
        if any(strcmp(tag,{'EXP_mono_tau_popup','EXP_ana_tau_popup'}))
          ResLibCal_UpdateDTau(h);
        elseif any(strcmp(tag,{'EXP_efixed','EXP_Kfixed','EXP_Lfixed'}))
          ResLibCal_UpdateEKLfixed(h);
        end
      elseif strcmp(get(h,'Type'), 'uimenu')
        ResLibCal(get(h,'Tag'));
      elseif any(strcmp(tag,{'EXP_mono_d','EXP_ana_d'}))
        ResLibCal_UpdateTauPopup;
      end
      varargin(2)=[];
      % update computation and plots
      feval(mfilename, 'update');
    otherwise
      % open file name or list of parameters given as 'VAR=VAL; ...'  
      if numel(varargin) > 1 && isstruct(varargin{2})
        out = varargin{2};
        varargin(2)=[];
      end
      out = ResLibCal_Open(action, out); % update 'out/EXP' from file
      ResLibCal_EXP2fig(out);                        % put it into the main GUI
      if numel(varargin) == 0
        out = ResLibCal_Compute(out);                    % compute the resolution
        ResLibCal_UpdateViews(out); % update views when they exist
      end
      
    end % switch (action)
    varargin(1) = [];
    % end if varargin is char
  elseif isstruct(varargin{1})
    % read an out or EXP structure
    EXP = varargin{1};
    if isfield(EXP,'EXP'),
      out = EXP; EXP=out.EXP;
    end
    varargin(1) = [];
  elseif numel(varargin) >= 4 && isnumeric(varargin{1}) && isnumeric(varargin{2}) ...
    && isnumeric(varargin{3}) && isnumeric(varargin{4}) 
    % read HKLE coordinates and compute resolution there
    % get current config
    fig = ResLibCal_fig;
    if ~isempty(fig)
      EXP = ResLibCal_fig2EXP(fig);
      if isfield(EXP,'EXP'),
        out = EXP; EXP=out.EXP;
      end
    end
    if isempty(out)
      filename = fullfile(prefdir, 'ResLibCal.ini');
      out = ResLibCal_Open(filename); % open the 'ResLibCal.ini' file (last saved configuration)
      if isfield(out,'EXP'),
        out = EXP; EXP=out.EXP;
      else
        EXP=out;
      end
    end
    if isempty(EXP) || ~isstruct(EXP), return; end

    if ~isempty(varargin) && isnumeric(varargin{1}) 
      if ~isempty(varargin{1}),
        EXP.QH = varargin{1};
        set(ResLibCal_fig('EXP_QH'),'String', mat2str(EXP.QH));
      end
      varargin(1)=[];
    end
    if ~isempty(varargin) && isnumeric(varargin{1}) 
      if ~isempty(varargin{1}),
        EXP.QK = varargin{1};
        set(ResLibCal_fig('EXP_QK'),'String', mat2str(EXP.QK));
      end
      varargin(1)=[];
    end
    if ~isempty(varargin) && isnumeric(varargin{1}) 
      if ~isempty(varargin{1}),
        EXP.QL = varargin{1};
        set(ResLibCal_fig('EXP_QL'),'String', mat2str(EXP.QL));
      end
      varargin(1)=[];
    end
    if ~isempty(varargin) && isnumeric(varargin{1}) 
      if ~isempty(varargin{1}),
        EXP.W = varargin{1};
        set(ResLibCal_fig('EXP_W'),'String', mat2str(EXP.W));
      end
      varargin(1)=[];
    end
    if isempty(varargin)
      out = ResLibCal_Compute(EXP);
    end
    if ~isempty(fig), ResLibCal_UpdateViews(out); end % when they exist
  end % if type(varargin)
  
  
end % end while nargin > 0
% end ResLibCal main

% ==============================================================================
% most functions in 'private'
%

% ==============================================================================
function filename = ResLibCal_Save
  filename = ResLibCal_Saveas(fullfile(prefdir, 'ResLibCal.ini'));

% ==============================================================================
function out = ResLibCal_UpdateViews(out)
% ResLibCal_UpdateViews: update all views (only when already visible)
%
  if nargin == 0, out = ''; end
  if ~isstruct(out), out = ResLibCal_Compute; end
  fig = ResLibCal_fig;
  if ~isempty(fig) 
    if strcmp(get(ResLibCal_fig('View_AutoUpdate'), 'Checked'), 'on')
      out = ResLibCal_UpdateResolution1(out); % TAS geometry
      out = ResLibCal_UpdateResolution2(out); % 2D, also shows matrix
      out = ResLibCal_UpdateResolution3(out); % 3D
    end
    ResLibCal_MethodEnableControls(out);    % enable/disable widgtes depending on the chosen method
  end
  % if no view exists, send result to the console 
  % here unactivated in case we use it as a model for e.g. fitting
  if isempty(fig) || isempty([ findobj(0, 'Tag','ResLibCal_View2') ...
    findobj(0, 'Tag','ResLibCal_View3') ])
		% display result in the console
		rlu = get(ResLibCal_fig('View_ResolutionRLU'), 'Checked');
		if ~strcmp(rlu, 'on'), mode=''; else mode='rlu'; end
		[res, inst] = ResLibCal_FormatString(out, mode);
		disp(char(res));
		disp(char(inst));
  end

% ==============================================================================
function out = ResLibCal_ViewResolution(out, dim)
% ResLibCal_ViewResolution: open the Resolution 2D/3D plot view
%
  if nargin == 0, out = ''; end
  if ~isstruct(out), out = ResLibCal_Compute; end
  h = findobj(0, 'Tag',[ 'ResLibCal_View' num2str(dim)]);
  if isempty(h)
    if dim~=1, name=sprintf('(%iD)', dim); else name='Matrix'; end
    h = figure('Name',[ 'ResLibCal: View Resolution ' name ], ...
               'Tag', [ 'ResLibCal_View' num2str(dim)], 'ToolBar','figure');
    p = get(h, 'Position'); p(3:4) = [ 640 480 ]; set(h, 'Position',p);
  else
    figure(h);
  end

% ==============================================================================
function out = ResLibCal_UpdateResolution1(out)
% ResLibCal_UpdateResolution1: update the TAS geometry view
%
  if nargin == 0, out = ''; end
  if ~isstruct(out), out = ResLibCal_Compute; end
  h = findobj(0, 'Tag','ResLibCal_View1');
  if isempty(h), return; end
  set(0,'CurrentFigure', h);
  set(h, 'Name','ResLibCal: View TAS geometry');

  % update/show the TAS geometry
  out = ResLibCal_TASview(out);

% ==============================================================================
function out = ResLibCal_UpdateResolution2(out)
% ResLibCal_UpdateResolution2: update the 2D view
%
  if nargin == 0, out = ''; end
  if ~isstruct(out), out = ResLibCal_Compute; end
  h = findobj(0, 'Tag','ResLibCal_View2');
  if isempty(h), return; end
  set(0,'CurrentFigure', h);

  % update/show the resolution projections
  rlu = get(ResLibCal_fig('View_ResolutionRLU'), 'Checked');
  qz  = get(ResLibCal_fig('View_ResolutionXYZ'), 'Checked');
  if strcmp(rlu, 'on'), rlu='rlu'; end
  if strcmp(qz, 'on'),  qz='qz'; end
  out = ResLibCal_Plot2D(out, [ rlu ' ' qz ]);

function out = ResLibCal_UpdateResolution3(out)
% ResLibCal_UpdateResolution3: update the 3D view
%
  if nargin == 0, out = ''; end
  if ~isstruct(out), out = ResLibCal_Compute; end
  h = findobj(0, 'Tag','ResLibCal_View3');
  if isempty(h), return; end
  set(0,'CurrentFigure', h);

  % update/show the resolution projections
  rlu = get(ResLibCal_fig('View_ResolutionRLU'), 'Checked');
  qz  = get(ResLibCal_fig('View_ResolutionXYZ'), 'Checked');
  if strcmp(rlu, 'on'), rlu='rlu'; end
  if strcmp(qz, 'on'),  qz='qz'; end
  out = ResLibCal_Plot3D(out, [ rlu ' ' qz ]);

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
      if ~isempty(strfind(lower(out.EXP.method), 'popovici')) % this is Popovici method
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
  
  
