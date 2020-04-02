function out = load_NeXus(in)
% LOAD_NEXUS format an HDF file with NeXus structure (signal, axes, ...)
%   See https://www.nexusformat.org/
%
% input:
%   in:  initial single HDF/NeXus data set loaded as a raw iData
% returns:
%   out: NXdata and NXdetector blocks
%
% called by: openhdf

% The elemnets we wish to identify and expose.
%   NXdetector
%   NXdata
%   NXsubentry
%   NXentry
%   NXinstrument
%   NXsample
%   NXprocess
%   NXuser
%
% Attributes can be as follows:
%   NXdata@signal= '<signal>'             name of field with signal (string)
%   NXdata@axes  = '<axis1>,<axis2>,...'  name of axes, in order (string, cellstr)
%                                         '.' indicates unspecified axis.
%
%   NXdata.<signal>@signal=1          indicates this is the signal
%   NXdata.<signal>@axes='<axis1>,<axis2>,...' name of axes for the signal
%   NXdata.<axis1>@axis=<rank>          indicates this is an axis of given rank

out     = [];

% We search for the 'NX_class' items. These identify the entries.
% These are defined as Attributes.
[NX_class, NX_path] = load_NeXus_class(in);
if isempty(NX_class), return; end % Not a NeXus file

% location of NXdata stuff
NXdata_path = NX_path(strcmp(NX_class, 'NXdata')  | strcmp(NX_class, 'NXdetector'));
if isempty(NXdata_path), return; end
% location of NXentry blocks
NXentry_path= NX_path(strcmp(NX_class, 'NXentry') | strcmp(NX_class, 'NXsubentry')); 

% get the Attributes to NXdata groups
NXdata_attr = fileattrib(in, NXdata_path);
if ~iscell(NXdata_path), NXdata_path = { NXdata_path }; end
if ~iscell(NXdata_attr), NXdata_attr = { NXdata_attr }; end

% We search for 'signal' Attributes (exact case). This defines the Signal.
NXdata_path_signal = [];
for index=1:numel(NXdata_path)
  if ~isstruct(NXdata_attr{index}) || isempty(NXdata_attr{index}), continue; end
  % we use the overloaded struct.findfield shipped with iFit to search for 'signal'
  if  isfield(          NXdata_attr{index}, 'signal') ...
  || ~isempty(findfield(NXdata_attr{index}, 'signal','exact case'))
    % this is an NXdata with 'signal' Attribute.
    NXdata_path_signal(end+1) = index; % NXdata index which has a signal inside.
  end
end

if in.verbose > 1
  disp([ mfilename ': found ' num2str(numel(NXdata_path_signal)) ...
    ' data sets with Signal in object ' in.Tag ]);
end
if isempty(NXdata_path_signal), return; end % No Signal in any NXdata/NXdetector

% we create as many output objects as NXdata and NXdetector entries
all_signal_names = {};
for index=1:numel(NXdata_path_signal)
  % all assignments are sent to the end of the loop to avoid intermediate checks.
  this_NXdata = NXdata_path{NXdata_path_signal(index)}; % NXdata path containing a Signal

  [this,this_signal_name] = load_NeXus_single(in, ...
    this_NXdata, NXentry_path, all_signal_names);
  % add a new data set when Signal is not already found.
  if ~isempty(this)
    out = [ out this ];
    all_signal_names{end+1} = this_signal_name;
  end

end

% ------------------------------------------------------------------------------
function [this, this_signal_name] = load_NeXus_single(in, ...
  this_NXdata, NXentry_path, all_signal_names)
% LOAD_NEXUS_SINGLE Analyse the attributes and determine the NXentry above
%   the given NXdata, as well as the inner Signal, Axes, and other symbols.

  this = copyobj(in); % we work on a copy
  if isempty(this_NXdata) || ~ischar(this_NXdata), return; end

  % determine the NXentry we are in (when multiple) ============================
  % The NXentry_path corresponds with the beginning of the this_NXdata
  % where the signal resides.
  this_NXentry = []; NXentry_remove = {};
  for index_entry=1:numel(NXentry_path)
    % store NXentry in which the Signal is located.
    if strncmp(NXentry_path{index_entry}, this_NXdata, numel(NXentry_path{index_entry}))
      this_NXentry = NXentry_path{index_entry};
    else
      NXentry_remove{end+1} = NXentry_path{index_entry};
    end
  end
  % Then remove other NXentry groups and set an alias to the current NXentry
  if ~isempty(this_NXentry) && ~isempty(NXentry_remove)
    set(this, NXentry_remove, []); % follow links to empty the targets
  end
  % we use 'NX_' preneding. Using bare 'NX' may clash with other entries in the object
  % so that e.g. 'NXdata' class may actually be resolved as a link in the object
  if ~isempty(this_NXentry), setalias(this, 'NX_entry', this_NXentry); end
  setalias(this, 'NX_data',  this_NXdata);
  axescheck(this, 'nocheck'); % disable auto axes check
  if in.verbose > 1
    if ~isempty(this_NXentry), disp([ mfilename ': NXentry ' this_NXentry ]); end
    disp([ mfilename ':   NXdata '  this_NXdata  ]);
  end

  % get the NX_class in this lightweight object
  [NX_class_entry, NX_path_entry] = load_NeXus_class(this);

  % remove all other NXdata items ==============================================
  % this allows to properly search within a single NXdata, without any chance of
  % interference with other data blocks.
  index_all = strcmp('NXdata', NX_class_entry) | strcmp('NXdetector', NX_class_entry);
  index_this= find(strncmp(this_NXdata, NX_path_entry, numel(this_NXdata)));
  index_all(index_this) = 0; % make sure we do not remove the active NXdata block
  % There is an issue here: The actual signal can be removed when set as an inner NXdata.
  set(this, NX_path_entry(find(index_all)), []); % follow links to empty targets
  axescheck(this, 'nocheck'); % disable auto axes check

  % determine the Signal =======================================================
  % Attributes can be as follows:
  %   NXdata@signal= '<signal>'         name of field with signal (string)
  %   NXdata.<signal>@signal=1          indicates this is the signal
  
  this_NXdata_attr = fileattrib(this, this_NXdata);
  this_signal_name = [];

  % is 'signal' an attribute of the NXdata ?
  if isfield(this_NXdata_attr, 'signal')
    %   NXdata@signal= '<signal>'         name of field with signal (string)
    this_signal_name = [ this_NXdata '.' this_NXdata_attr.signal ];
  else
    %   NXdata.<signal>@signal=1          indicates this is the signal (old style)
    % the attribute is attached to the Signal itself. Search it...
    attr_signal_path = findfield(this, 'signal','exact case');
    % search for a '.signal' that is an Attribute.
    isattribute      = cellfun(@numel, strfind(attr_signal_path, '.Attributes'));
    attr_signal_path = attr_signal_path(find(isattribute));
    % get the 'signal' attribute value: Is it a boolean (1) ?
    if ~isempty(attr_signal_path)
      if iscell(attr_signal_path), attr_signal_path = attr_signal_path{1}; end
      attr_signal_value = get(this, attr_signal_path);
      if isnumeric(attr_signal_value) && isscalar(attr_signal_value) ...
        && attr_signal_value
        % clean Attribute path. Remove 'Attributes' and get the name before '.signal'.
        attr_signal_path = strrep(attr_signal_path, '.Attributes', '');
        index = find(attr_signal_path == '.', 1, 'last');
        this_signal_name = attr_signal_path(1:(index-1));
      end
    end
  end
  % was this Signal found before ?
  if ~isempty(strcmp(this_signal_name, all_signal_names))
    this = []; return
  end
  
  % check if the Signal is indeed numeric
  if ~isempty(this_signal_name) && ~isempty(findfield(in, this_signal_name, 'numeric'))
    setalias(this, 'Signal',  this_signal_name);
    axescheck(this, 'nocheck'); % disable auto axes check
    if in.verbose > 1
      disp([ mfilename ':     Signal ' this_signal_name  ]);
    end
    % assign a Label to the Signal, using the Attributes.
    attr_signal_value = fileattrib(this, 'Signal');
    if ~isempty(attr_signal_value)
      if isfield(attr_signal_value, 'long_name')
        attr_signal_value = attr_signal_value.long_name;
      else
        attr_signal_value = class2str(' ', attr_signal_value, 'no comment short');
      end
      label(this, 'Signal', attr_signal_value);
    end
  else
    if in.verbose > 1
      disp([ mfilename ':     WARNING: no numeric Signal found in this NXdata.' ]);
    end
    return
  end

  % define 'NX' aliases in this Data block =====================================
  % NXinstrument NXprocess NXsample may be used from top-level. Expose them.
  for l={'instrument','sample','process','user'}
    NXblock = NX_path_entry(strcmp(NX_class_entry, [ 'NX' l{1} ]));
    if ~isempty(NXblock) && ~isfield(this, l{1})
      % Alias: e.g. instrument -> NXinstrument block
      setalias(this, [ 'NX_' l{1} ], NXblock{1});
      if in.verbose > 1
        disp([ mfilename ': NX' l{1} ' ' NXblock{1} ]);
      end
    end
  end
  axescheck(this, 'nocheck'); % disable auto axes check

  % determine the axes =========================================================
  % Attributes can be as follows:
  %   NXdata@axes  = '<axis1>,<axis2>,...'  name of axes, in order (string, cellstr)
  %                                         '.' indicates unspecified axis.
  %   NXdata.<signal>@axes='<axis1>,<axis2>,...' name of axes for the signal
  %   NXdata.<axis1>@axis =<rank>          indicates this is an axis of given rank
  this_axes_names = {};
  add_NXdata_path = true; % flag to prepend NXdata path 
  if isfield(this_NXdata_attr, 'axes')
    %   NXdata@axes  = '<axis1>,<axis2>,...'  name of axes, in order (string, cellstr)
    this_axes_names = this_NXdata_attr.axes;
  elseif ~isempty(findfield(this, 'axes','exact case'))
    %   NXdata.<signal>@axes='<axis1>,<axis2>,...' name of axes for the signal
    attr_axes_path = findfield(this, 'axes','exact case');
    % search for a '.axes' that is an Attribute.
    isattribute    = cellfun(@numel, strfind(attr_axes_path, '.Attributes'));
    attr_axes_path = attr_axes_path(find(isattribute));
    % get the 'axes' attribute value: is it a string/cellstr ?
    if ~isempty(attr_axes_path)
      if iscell(attr_axes_path), attr_axes_path = attr_axes_path{1}; end
      this_axes_names= get(this, attr_axes_path);  % a string/cellstr with axes names
    end
  elseif ~isempty(findfield(this, 'axis','exact case'))
    %   NXdata.<axis1>@axis =<rank>          indicates this is an axis of given rank
    attr_axes_path = findfield(this, 'axis','exact case');
    for index=1:numel(attr_axes_path)
      attr_axes_value= get(this, attr_axes_path{index}); % a numeric/scalar
      if isnumeric(attr_axes_value) && isscalar(attr_axes_value)
        % clean Attribute path. Remove 'Attributes' and get the name before '.axis'.
        attr_axis_path = strrep(attr_axes_path{index}, '.Attributes', '');
        index = find(attr_axis_path == '.', 1, 'last');
        this_axes_names{attr_axes_value} = attr_axis_path(1:(index-1)); % full path
        add_NXdata_path = false;
      end
    end
  end
  if ~isempty(this_axes_names) && add_NXdata_path
    if ischar(this_axes_names)
      this_axes_names = textscan(this_axes_names, '%s','Delimiter',','); % split in tokens
      this_axes_names = this_axes_names{1};
    end
    this_axes_names = strcat(this_NXdata, '.', this_axes_names); % NXdata.<axes>
  end
  if ischar(this_axes_names), this_axes_names = cellstr(this_axes_names); end
  % assign axes when these are numeric
  for index=1:numel(this_axes_names)
    this_name = this_axes_names{index};
    if isempty(this_name), continue; end
    if this_name(end) ~= '.' && ~isempty(findfield(in, this_name, 'numeric'))
      setaxis(this, index, this_name);
    end
    % assign a Label to the Axis, using the Attributes.
    attr_axis_value = fileattrib(this, index);
    if ~isempty(attr_axis_value)
      if isfield(attr_axis_value, 'long_name')
        attr_axis_value = attr_axis_value.long_name;
      else
        attr_axis_value = class2str(' ', attr_axis_value, 'no comment short');
      end
      label(this, index, attr_axis_value);
    end
  end
  if in.verbose > 1
    for index=1:numel(this_axes_names)
      disp([ mfilename ':     Axis{' num2str(index) '} ' this_axes_names{index} ]);
    end
  end
  

% ------------------------------------------------------------------------------
function [NX_class, NX_path] = load_NeXus_class(in)
% LOAD_NEXUS_CLASS Determine the NeXus classes and their path as data sets
  NX_path  = [];
  NX_class = [ findfield(in, 'NX_class', 'exact case') ...
               findfield(in, 'class', 'exact case') ];
  if isempty(NX_class), return; end % Not a NeXus file

  % Then we remove the '.Attributes' and '.NX_class' tokens to get the true path.
  NX_path = strrep(NX_class, '.Attributes','');
  NX_path = strrep(NX_path,  '.NX_class',  ''); % path to NX entries
  NX_path = strrep(NX_path,  '.class',  '');    % path to NX entries (old style)
  NX_class= get(in, NX_class);                  % type of NX entries (NX_class)
