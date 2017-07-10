function out = ResLibCal_UpdateViews(out, modev)
% ResLibCal_UpdateViews: update all views (only when already visible)
% modev can be: 'force' (update all views) or 'stdout'
%
  if nargin == 0, out = ''; end
  if nargin < 2, modev=''; end
  if ~isstruct(out), out = ResLibCal_Compute; end
  fig = ResLibCal_fig;
  if ~isempty(fig) || strcmp(modev, 'force')
    if strcmp(get(ResLibCal_fig('View_AutoUpdate'), 'Checked'), 'on') || strcmp(modev, 'force')
      t=clock;
      out = ResLibCal_UpdateResolution1(out); % TAS geometry
      out = ResLibCal_UpdateResolution2(out); % 2D, also shows matrix
      out = ResLibCal_UpdateResolution3(out); % 3D
      if ~strcmp(modev, 'force') && etime(clock, t) > 5
        disp([ mfilename ': the time required to update all plots gets long.' ])
        disp('INFO          Setting View/AutoUpdate to off.')
        set(ResLibCal_fig('View_AutoUpdate'), 'Checked', 'off');
      end
    end
    ResLibCal_MethodEnableControls(out);    % enable/disable widgtes depending on the chosen method
  end
  % if no view exists, send result to the console 
  % here unactivated in case we use it as a model for e.g. fitting
  if isempty(fig) || strcmp(modev, 'stdout') ...
  || isempty([ findobj(0, 'Tag','ResLibCal_View2') findobj(0, 'Tag','ResLibCal_View3') ])
		% display result in the console
		rlu = get(ResLibCal_fig('View_ResolutionRLU'), 'Checked');    % [a* b*  c* ]
		spec= get(ResLibCal_fig('View_ResolutionSPEC'),'Checked');    % [Ql Qt  Qv ]
		abc = get(ResLibCal_fig('View_ResolutionABC'), 'Checked');    % [A  B   C  ]
		lat = get(ResLibCal_fig('View_ResolutionLattice'), 'Checked');% [a* b'* c'*]
		modev='abc'; % default
		if     strcmp(rlu, 'on') modev='rlu'; 
		elseif strcmp(spec,'on') modev='spec'; 
		elseif strcmp(abc, 'on') modev='abc';
		elseif strcmp(lat, 'on') modev='lattice'; end
		[res, inst] = ResLibCal_FormatString(out, modev);
		disp(char(res));
		disp(char(inst));
  end
