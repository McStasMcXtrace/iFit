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
%   open, save, saveas, export, exit, reset, print, create,
%   compute, update (=compute+show), view2, view3, view_tas
%   default, quit, close
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
% Version: $Revision$ $Date$

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

ResLibCal_version = [ mfilename ' $Revision$ ($Date$)' ];

% menu actions:
if ~isempty(varargin)
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
    case {'file_open','open'}
      EXP = ResLibCal_fig2EXP(get(0,'CurrentFigure'));
      if isfield(EXP,'EXP'), EXP = EXP.EXP; end
      if length(varargin) < 1, varargin{2} = ''; end
      if length(varargin) < 2, varargin{3} = EXP; end
      out = ResLibCal_Open(varargin{2:end});  % (filename, EXP)
      out = ResLibCal('update');
    case 'file_open_instr'
      p = fileparts(which(mfilename));
      ResLibCal('open',[ p filesep 'instruments' ]);
    case {'file_reset','reset'}
      filename = fullfile(prefdir, 'ResLibCal.ini');
      if ~exist(filename, 'file')
        source = 'the factory defaults.';
      else
        source = [ 'the file ' filename ' (delete it to return to factory defaults)' ];
      end
      options.Default     = 'No';
      options.Interpreter = 'tex';
      ButtonName = questdlg({ ...
        '{\fontsize{14}{\color{blue}Reload default configuration ?}}', ...
        'This will reset all fields of the ResLibCal window from ', ...
        source, ...
        '{\bf Reset now ?}'}, 'ResLibCal: Reset ?', 'Yes', 'No', options);
      if strcmp(ButtonName, 'Yes')
        if exist(filename,'file')
          feval(mfilename, 'create');
        else
          fig = findall(0, 'Tag','ResLibCal');
          delete(fig);
          openfig('ResLibCal');
        end
        % make sure the QH,QK,QL,EN are all scalars
        EXP = ResLibCal_fig2EXP;
        EXP.QH=EXP.QH(1);
        EXP.QK=EXP.QK(1);
        EXP.QL=EXP.QL(1);
        EXP.W =EXP.W (1);
        ResLibCal_EXP2fig(EXP);
        out = ResLibCal('update');
      end
    case {'file_save','save'}
      % save configuration so that it is re-opened at next re-start
      ResLibCal_Save; % (filename=prefdir)
    case {'file_saveas','saveas'}
      % save configuration
      ResLibCal_Saveas(varargin{2:end}); % (filename, EXP)
    case {'file_print','print'}
      fig = findall(0, 'Tag','ResLibCal');
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
      fig = findall(0, 'Tag','ResLibCal');
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
      if length(varargin) > 1, v=varargin{2:end}; else v=[]; end
      out = ResLibCal_ViewResolution(v,2);  % open/raise View Res2
      out = ResLibCal_UpdateViews(out);
    case {'view_resolution3','view3'}
      if length(varargin) > 1, v=varargin{2:end}; else v=[]; end
      out = ResLibCal_ViewResolution(v,3);  % open/raise View Res2
      out = ResLibCal_UpdateViews(out);
    case {'view_tas','geometry','tas'}
      if length(varargin) > 1, v=varargin{2:end}; else v=[]; end
      out = ResLibCal_ViewResolution(v,1);  % open/raise View TAS
      out = ResLibCal_UpdateViews(out);
    case 'help_content'
      link = fullfile(fileparts(which(mfilename)), 'doc', [ mfilename '.html' ]);
      disp([ mfilename ': opening help from ' link ])
      web(link);
    case 'help_about'
      % get the ILL logo from object
      fig = findall(0, 'Tag','ResLibCal');
      if isempty(fig), return; end
      cdata = get(findall(fig, 'Tag', 'ILL_logo'), 'CData');
      message = {...
        [ '{\fontsize{14}{\color{blue}' ResLibCal_version '} EUPL license} ' ], ...
        'ResLibCal is a graphical user interface to compute and display' , ...
        'the {\bf triple-axis resolution function} obtained from e.g. Cooper-Nathans and Popovici analytical approximations. The GUI allows to select among a set of computation kernels.' , ...
        'This application was written by E. Farhi {\copyright}ILL/DS/CS <farhi@ill.eu> using' , ...
        '\bullet ResLib 3.4 (A. Zheludev)' , ...
        '\bullet ResCal5 (A. Tennant and D. Mc Morrow)' , ...
        '\bullet Res3ax (J. Ollivier)' , ...
        '\bullet Rescal/AFILL (Hargreave, Hullah) and vTAS view (A. Bouvet/A. Filhol)' };
      CreateMode.WindowStyle = 'modal';
      CreateMode.Interpreter='tex';
      msgbox(message, ...
        'About: ResLibCal', ...
        'custom',cdata,jet,CreateMode);
    case 'view_update'
      out = ResLibCal_Compute(varargin{2:end}); % arg can be an EXP
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
    % other actions (not menu items) -------------------------------------------
    case {'create','default'}
      fig = findall(0, 'Tag','ResLibCal');
      if isempty(fig) || ~ishandle(fig)
        disp([ 'Welcome to ' ResLibCal_version ]);
        openfig('ResLibCal'); % open the main ResLibCal figure.
        if strcmp(action, 'create')
          filename = fullfile(prefdir, 'ResLibCal.ini');
          out = ResLibCal_Open(filename); % open the 'ResLibCal.ini' file (last saved configuration)
          out = ResLibCal('compute');
        end
      elseif length(fig) > 1
        delete(fig(2:end)); % remove duplicated windows
      end
      out = ResLibCal_fig2EXP(fig);
    case 'update'
      % update all opened views with new computation (widget update)
      fig = findall(0, 'Tag','ResLibCal');
      
      out = feval(mfilename, 'compute');
      
      if ~isempty(fig) && strcmp(get(findobj(fig, 'Tag','View_AutoUpdate'), 'Checked'), 'on')
        out = ResLibCal_UpdateViews(out);
      elseif isempty(fig)
        out = ResLibCal_UpdateViews(out);
      end
    case {'compute','resolution'}
      % only compute. No output except in varargout
      fig = findall(0, 'Tag','ResLibCal');
      % if no interface exists, load the last saved configuration before computing
      if isempty(fig)
        filename = fullfile(prefdir, 'ResLibCal.ini');
        out = ResLibCal_Open(filename); % open the 'ResLibCal.ini' file (last saved configuration)
        ResLibCal_UpdateViews; % when they exist
      else
        out = ResLibCal_Compute(varargin{2:end}); % arg can be an EXP
      end
    case {'view_close','close'}
      delete(findall(0, 'Tag', 'ResLibCal_View1'));
      delete(findall(0, 'Tag', 'ResLibCal_View2'));
      delete(findall(0, 'Tag', 'ResLibCal_View3'));
    case 'update_d_tau'
      % update d-spacing from a popup item
      ResLibCal_UpdateDTau(varargin{2:end});      % arg is popup handle
    case 'update_d_tau_popup'
      % update the mono-ana popup from the DM/DA when value is close
      ResLibCal_UpdateTauPopup;
    case 'update_ekl'
      % update E, K, lambda
      ResLibCal_UpdateEKLfixed(varargin{2:end});  % arg is edit box handle
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
      % update computation and plots
      feval(mfilename, 'update');
    otherwise
      try
        ResLibCal;
        out = ResLibCal_Open(action, varargin{2:end});
        ResLibCal_EXP2fig(out);
        out = ResLibCal('compute');
        ResLibCal_UpdateViews; % when they exist
      catch
        disp([ mfilename ': Unknown action ' action ]);
        out = [];
      end
    end
    % end if varargin is char

  elseif nargin >=1
    if isstruct(varargin{1})  % an EXP structure ?
      EXP = varargin{1};
      varargin(1)=[];
      fig = [];
    else
      [EXP, fig] = ResLibCal_fig2EXP;
    end
    if isfield(EXP,'EXP'),
      out = EXP; EXP=out.EXP;
    end
    if isempty(EXP) || ~isstruct(EXP), return; end

    if ~isempty(varargin) && isnumeric(varargin{1}) 
      if ~isempty(varargin{1}),
        EXP.QH = varargin{1};
        set(findall(fig,'Tag','EXP_QH'),'String', mat2str(EXP.QH));
      end
      varargin(1)=[];
    end
    if ~isempty(varargin) && isnumeric(varargin{1}) 
      if ~isempty(varargin{1}),
        EXP.QK = varargin{1};
        set(findall(fig,'Tag','EXP_QK'),'String', mat2str(EXP.QK));
      end
      varargin(1)=[];
    end
    if ~isempty(varargin) && isnumeric(varargin{1}) 
      if ~isempty(varargin{1}),
        EXP.QL = varargin{1};
        set(findall(fig,'Tag','EXP_QL'),'String', mat2str(EXP.QL));
      end
      varargin(1)=[];
    end
    if ~isempty(varargin) && isnumeric(varargin{1}) 
      if ~isempty(varargin{1}),
        EXP.W = varargin{1};
        set(findall(fig,'Tag','EXP_W'),'String', mat2str(EXP.W));
      end
      varargin(1)=[];
    end
    out = ResLibCal_Compute(EXP);
    if ~isempty(fig), ResLibCal_UpdateViews(out); end % when they exist
  end
  % end nargin > 0
else
  % nargin == 0

  % open or create the main GUI
  feval(mfilename, 'create'); % load last configuration
  out = ResLibCal_Compute;
  out = ResLibCal_ViewResolution(out,2);  % open/raise View Res2
  out = ResLibCal_UpdateViews(out); % when they exist
end
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
  out = ResLibCal_UpdateResolution1(out); % TAS geometry
  out = ResLibCal_UpdateResolution2(out); % 2D, also shows matrix
  out = ResLibCal_UpdateResolution3(out); % 3D
  % if no view exists, send result to the console 
  % here unactivated in case we use it as a model for e.g. fitting
  if isempty([ findall(0, 'Tag','ResLibCal_View2') ...
    findall(0, 'Tag','ResLibCal_View3') ])
		% display result in the console
		[res, inst] = ResLibCal_FormatString(out);
		disp(char(res));
		disp(char(inst));
  end

% ==============================================================================
function out = ResLibCal_ViewResolution(out, dim)
% ResLibCal_ViewResolution: open the Resolution 2D/3D plot view
%
  if nargin == 0, out = ''; end
  if ~isstruct(out), out = ResLibCal_Compute; end
  h = findall(0, 'Tag',[ 'ResLibCal_View' num2str(dim)]);
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
  h = findall(0, 'Tag','ResLibCal_View1');
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
  h = findall(0, 'Tag','ResLibCal_View2');
  if isempty(h), return; end
  set(0,'CurrentFigure', h);

  % update/show the resolution projections
  rlu = get(findobj(out.handle,'Tag','View_ResolutionRLU'), 'Checked');
  qz  = get(findobj(out.handle,'Tag','View_ResolutionXYZ'), 'Checked');
  if strcmp(rlu, 'on'), rlu='rlu'; end
  if strcmp(qz, 'on'),  qz='qz'; end
  out = rc_projs(out, [ rlu ' ' qz ]);

function out = ResLibCal_UpdateResolution3(out)
% ResLibCal_UpdateResolution3: update the 3D view
%
  if nargin == 0, out = ''; end
  if ~isstruct(out), out = ResLibCal_Compute; end
  h = findall(0, 'Tag','ResLibCal_View3');
  if isempty(h), return; end
  set(0,'CurrentFigure', h);

  % update/show the resolution projections
  % update/show the resolution projections
  rlu = get(findobj(out.handle,'Tag','View_ResolutionRLU'), 'Checked');
  qz  = get(findobj(out.handle,'Tag','View_ResolutionXYZ'), 'Checked');
  if strcmp(rlu, 'on'), rlu='rlu'; end
  if strcmp(qz, 'on'),  qz='qz'; end
  out = ResPlot3D(out, [ rlu ' ' qz ]);

function ResLibCal_UpdateTauPopup
% update the popup menu from the editable mono/ana value when d is close
%
  [EXP,fig] = ResLibCal_fig2EXP;
  if isempty(fig), return; end
  popup = findobj(fig,'Tag','EXP_mono_tau_popup');
  label = GetTau(EXP.mono.tau, 'getlabel');
  index = find(strncmpi(get(popup,'String'), label, length(label)));
  if ~isempty(index) && ~isempty(label), set(popup, 'value', index(1)); end

  popup = findobj(fig,'Tag','EXP_ana_tau_popup');
  label = GetTau(EXP.ana.tau, 'getlabel');
  index = find(strncmpi(get(popup,'String'), label, length(label)));
  if ~isempty(index) && ~isempty(label), set(popup, 'value', index(1)); end
