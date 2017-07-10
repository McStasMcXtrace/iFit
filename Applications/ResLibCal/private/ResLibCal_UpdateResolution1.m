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

