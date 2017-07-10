function out = ResLibCal_UpdateResolution3(out)
% ResLibCal_UpdateResolution3: update the 3D view
%
  if nargin == 0, out = ''; end
  if ~isstruct(out), out = ResLibCal_Compute; end
  h = findobj(0, 'Tag','ResLibCal_View3');
  if isempty(h), return; end
  set(0,'CurrentFigure', h);

  % update/show the resolution projections
  rlu = get(ResLibCal_fig('View_ResolutionRLU'), 'Checked');    % [a* b*  c* ]
	spec= get(ResLibCal_fig('View_ResolutionSPEC'),'Checked');    % [Ql Qt  Qv ]
	abc = get(ResLibCal_fig('View_ResolutionABC'), 'Checked');    % [A  B   C  ]
	lat = get(ResLibCal_fig('View_ResolutionLattice'), 'Checked');% [a* b'* c'*]
	if     strcmp(rlu, 'on') modev='rlu'; 
	elseif strcmp(spec,'on') modev='spec'; 
	elseif strcmp(abc, 'on') modev='abc';
	elseif strcmp(lat, 'on') modev='lattice'; end
	
  qz  = get(ResLibCal_fig('View_ResolutionXYZ'), 'Checked');
  qyz = get(ResLibCal_fig('View_ResolutionHV'), 'Checked');
  MC  = get(ResLibCal_fig('View_Resolution_Cloud'), 'Checked');
  if strcmp(qz, 'on'),  qz='qz'; else  qz = 'en'; end
  if strcmp(qyz,'on'), qyz='xz'; else qyz = 'xy'; end
  if strcmp(MC, 'on'),  MC='cloud'; end
  out = ResLibCal_Plot3D(out, [ modev ' ' qyz ' ' qz ' ' MC ]);
