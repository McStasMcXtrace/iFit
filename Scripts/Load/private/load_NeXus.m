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

out = in;
verbose = in.verbose;

% We search for the 'NX_class' items. These identify the entries.
% These are defined as Attributes.
NX_class = findfield(in, 'NX_class');
if isempty(NX_class), return; end % Not a NeXus file

% Then we remove the '.Attributes' and '.NX_class' tokens to get the true path.
NX_path = strrep(NX_class, '.Attributes','');
NX_path = strrep(NX_path,  '.NX_class',  ''); % path to NX entries
NX_class= get(in, NX_class);                  % type of NX entries (NX_class)

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
NXdata_path_signal = {};
for index=1:numel(NXdata_path)
  if ~isstruct(NXdata_attr{index}) || isempty(NXdata_attr{index}), continue; end
  % we use the overloaded struct.findfield shipped with iFit
  if ~isempty(findfield(NXdata_attr{index}, 'signal','exact case'))
    % this is an NXdata with 'signal' Attribute.
    NXdata_path_signal{end+1} = NXdata_path{index};
  end
end

if verbose > 1
  disp([ mfilename ': found ' num2str(numel(NXdata_path_signal)) ...
    ' data sets with Signal in object ' in.Tag ]);
end
if isempty(NXdata_path_signal), return; end % No Signal in any NXdata/NXdetector

% we create as many output objects as NXdata and NXdetector entries
for index=1:numel(NXdata_path_signal)
  % we define a new object with Signal
  if numel(NXdata_path_signal) > 1, this = copyobj(in); else this = in; end
  attr = NXdata_path{index};
  setaxis(this, 0, NXdata_path_signal{index});

  % assign a Label to the Signal, using the Attributes.
  lab = strrep(class2str(NXdata_path{index}, 'no comment short'), 'struct_str.','');
  label(  this, 0, strrep(lab, sprintf('\n'), '') );

  % we extract the NXdata which holds the Signal
  [~, this_NXdata] = fileparts(NXdata_path_signal{index});
  if verbose > 1
    disp([ mfilename ': defining Signal ' NXdata_path_signal{index} ' in object ' in.Tag ]);
    disp([ mfilename ':   NXdata ' this_NXdata ' with Signal in object ' in.Tag ]);
  end

  % determine the NXentry we are in (when multiple).
  % The NXentry_path corresponds with the beginning of the NXdata_path{index}
  % where the signal resides.
  this_NXentry = []; NXentry_remove = {};
  for index_entry=1:numel(NXentry_path)
    % store NXentry in which the Signal is located.
    if strncmp(NXentry_path{index_entry}, NXdata_path{index}, ...
        numel(NXentry_path{index_entry}))
      this_NXentry = NXentry_path{index_entry};
    else
      NXentry_remove{end+1} = NXentry_path{index_entry};
    end
  end
  % Then remove other NXentry groups and set an alias to the current NXentry
  if ~isempty(this_NXentry) && ~isempty(NXentry_remove)
    set(this, NXentry_remove, []);
  end
  set(this, 'NXentry', this_NXentry);
  if verbose > 1
    disp([ mfilename ':   NXentry ' this_NXentry ' with NXdata/Signal in object ' in.Tag ]);
  end

  % get the NX_class in this lightweight object
  NX_class_entry = findfield(this, 'NX_class');
  % remove the '.Attributes' and '.NX_class' tokens to get the true path.
  NX_path_entry = strrep(NX_class_entry, '.Attributes','');
  NX_path_entry = strrep(NX_path_entry,  '.NX_class',  ''); % path to NX entries
  NX_class_entry= get(this, NX_class_entry);                % type of NX entries (NX_class)
  
  % define 'NX' aliases in this Data block
  % NXinstrument NXprocess NXsample
  for l={'Instrument','Sample','Process','User'}
    NXblock = NX_path_entry(strcmp(NX_class_entry, [ 'NX' lower(l{1}) ]));
    if ~isempty(NXblock) && ~isfield(this, l{1})
      % Alias: e.g. instrument -> NXinstrument block
      setalias(this, l{1}, NXblock{1});
      if verbose > 1
        disp([ mfilename ':   NX ' l{1} ' as ' NXblock{1} ' in object ' in.Tag ]);
      end
    end
  end

  % now set the axes...
  % according to nexusformat.org the axes are defined with either:
  %   

  
  out = [ out this ];

end

