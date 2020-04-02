function out = openhdf(filename, format)
%OPENHDF Open an HDF file, display it
%        and set the 'ans' variable to an iData object with its content
% 

if nargin < 2
  format = 'HDF';
end

if ~isa(filename,'estruct') && ~isa(filename, 'iData')
  out = estruct(iLoad(filename,format));
else
  out = filename;
end
clear filename;

if numel(out) > 1
  % handle input iData arrays
  in = out;
  out = []; % the number of elements may change, can not simply replace
  for index=1:numel(in)
    out = [ out ; feval(mfilename, in(index)) ];
    in(index) = estruct; % free memory
  end
  return
end

if ~isempty(findstr(out, 'NeXus')) || ~isempty(findfield(out, {'NX_class','class'}))
  % special stuff for NeXus files
  out1 = load_NeXus(out); % see private
  
  % call other specific importers
  if ismethod(out1, 'findfield') && ~isempty(findfield(out1, 'RITA_2'))
    out1 = load_psi_RITA(out1); % see private
  end

  if isa(out1, 'estruct') || isa(out1, 'iData')
    out = out1;
  end
  
end % if Nexus
  
if ~nargout
  figure; subplot(out);
  
  if ~isdeployed
    assignin('base','ans',out);
    ans = out
  end
end

% ------------------------------------------------------------------------------

