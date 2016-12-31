function EXP = ResLibCal_Open(filename, EXP, silent)
% EXP = ResLibCal_Open(filename): open an EXP/ResLib file and update main GUI
%
% Input:
% filename: file name (ResLibCal, ResCal or ILL TAS data file), 
%           or character string to evaluate (producing a structure or a vector)
%           or character string describing a structure
%           or structure with either EXP ResLib or ResCal fields
%
% Return:
%  EXP: configuration structure

% Calls: ResLibCal_RescalPar2EXP, ResLibCal_EXP2fig

  if nargin < 1, filename = ''; end
  if nargin < 2, EXP = []; end
  if nargin < 3, silent = 0; end
  if isempty(filename) || (ischar(filename) && isdir(filename))
    [filename, pathname] = uigetfile( ...
       {'*.m;*.ini',  'ResLibCal configuration M-file (*.m;*.ini)'; ...
        '*.cfg;*.par;*.res','ResCal/ResCal5 configuration (*.par;*.cfg;*.res)' ; ...
        '*.*',  'All Files, including ILL TAS Data and LLB spec (*.*)'}, ...
        'Open configuration as ResLibCal, ResCal, ILL TAS Data, ...', filename);
    if isempty(filename) || all(filename == 0), return; end
    filename = fullfile(pathname, filename);
  end

  if ischar(filename) && exist(filename,'file') % a file exists: read it
    % handle case of ResCal5 .par .cfg file (numerical vector)
    try
      content = load(filename); % usually produces a vector of Rescal parameters (42 or 27)
                                % can also be a .mat file for EXP/out
    catch
      % if not ResCal vector, read the file content (char)
      content = fileread(filename);
    end
  else
    content = filename;
  end

  if ischar(content) && ~isempty(content)  % content of a file, or string to evaluate
    try
      evalc(content);% this should make an 'EXP' or 'config' variable
      if ~isempty(config)
        EXP = config; % replace EXP full config (override any previous setting)
        content = ''; % success unactivates further interpretation of input
      end
    end
  end
  titl = 'ResLibCal configuration';
  % converted from a string or read from a file
  if ischar(content) || isstruct(content) || isnumeric(content) 
    % read content as a structure, ResCal par/cfg, ...
    if isfield(EXP,'EXP'), EXP=EXP.EXP; end
    [EXP,titl] = ResLibCal_RescalPar2EXP(content, EXP);
    % overload EXP with ResCal structure if EXP is incomplete and ResCal is
    % there
    if isfield(EXP,'ResCal') && (~isfield(EXP, 'mono') || ~isfield(EXP, 'sample') || ~isfield(EXP, 'ana'))
      EXP = ResLibCal_RescalPar2EXP(EXP.ResCal, EXP);
    end
  end

  % evaluate it to get 'EXP'
  if isstruct(EXP)
    Position = [];
    if exist('config','var') && isstruct(config) 
      if isfield(config, 'Title')
        titl = config.Title;
      end
      if isfield(config, 'Position')
        Position = config.Position;
      end
    end
    % send it to the figure
    try
      [fig, EXP] = ResLibCal_EXP2fig(EXP); % open figure if not yet done
      % set position and size from last save (if available)
      if ~isempty(Position) && isnumeric(Position) && ~isempty(fig)
        if numel(Position) == 4
          set(fig, 'Units','pixels','Position',Position);
        elseif numel(Position) == 2 % size only
          set(fig, 'Units','pixels');
          p0 = get(fig,'Position');
          Position(1:2) = p0(1:2);
          set(fig, 'Position',Position);
        end
      end
      % force full update of all fields
      ResLibCal('update_d_tau_popup');
      if ~silent
        if isstruct(filename)
          disp([ datestr(now) ': Loaded ' titl ' from ' ]);
          disp(filename)
        else
          disp([ datestr(now) ': Loaded ' titl ' from ' filename ]);
        end
      end
    catch
      warning([ datestr(now) ': Could not load ResLibCal configuration '  ]);
      disp(filename)
      rethrow(lasterror)
    end
  end
% end ResLibCal_Open
