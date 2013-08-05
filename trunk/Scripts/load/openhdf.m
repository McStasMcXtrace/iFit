function out = openhdf(filename, format)
%OPENHDF Open an HDF file, display it
%        and set the 'ans' variable to an iData object with its content

if nargin < 2
  format = 'HDF';
end

if ~isa(filename,'iData')
  out = iData(iLoad(filename,format));
else
  out = filename;
end

if length(out(:)) > 1
  % handle input iData arrays
  for index=1:length(out(:))
    out(index) = feval(mfilename, out(index));
  end
else
  % convert all links to valid links in iData objects
  
  
  % links from hdf5 have path with structure.
  if ~isempty(findstr(out, 'NeXus')) || ~isempty(findstr(out, 'NX_class'))
    % special studff for NeXus files
    
    % identify Signal, and search for its Attributes: axes, signal=1

    % identify root Attributes

    % search for a 'process' group, and build a CommandHistory from it

    % get title, instrument.name
    
    % call other specific importers
    % load_psi_RITA
  end % if Nexus
end

if ~nargout
  figure; subplot(out);
  
  if ~isdeployed
    assignin('base','ans',out);
    ans = out
  end
end

