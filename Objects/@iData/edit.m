function edit(a)
% edit(s) : edit iData object Signal/Monitor
%
%   @iData/edit function to view the Signal/Monitor of a data set. 
%     The data appears as a spreadsheet (requires Java to be enabled).
%     The Signal can not (yet) be modified.
%
% input:  s: object or array (iData)
% output: b: object or array (iData)
% ex:     b=edit(a);
%
% Version: $Date$
% See also iData, iData/uminus, iData/abs, iData/real, iData/imag, iData/uplus

% the Selection is stored in the figure UserData

  if numel(a) > 1
    parfor index=1:numel(a)
      edit(a(index));
    end
    return
  end

  if exist('uitable')
    if prod(size(a)) > 1e5
      iData_private_warning(mfilename, [ 'Object ' a.Tag ' is too large (numel=' num2str(prod(size(a))) ...
    '.\n\tYou should rebin with e.g. a=a(1:2:end, 1:2:end, ...).' ]);
    end
    
    Signal    = getaxis(a, 'Signal'); % Signal/Monitor
    
    % opens a Table with the Signal content there-in
    NL = sprintf('\n');
    f = figure('Name', [ 'Edit ' char(a) ]); % raw Signal, no monitor weightening
    p = get(f, 'Position');
    h = uitable('Data', Signal);
    
    % UserData.Monitor = real(get(a,'Monitor'));
    
    UserData.Selection = [];
    UserData.handle    = h;
    UserData.Name      = char(a);
    
    set(h, 'Position', [0,0,p(3),p(4)]); 
    try
      set(h, 'Tag',[ mfilename '_' a.Tag ]);
      set(h, 'Units','normalized');
      set(h, 'Position', [0,0,1,1]);
      set(f, 'UserData', UserData);  % contains the selection indices
      set(h, 'CellSelectionCallback', @iData_edit_CellSelectionCallback);
    end

%      'CellEditCallback', @iData_edit_CellEditCallback, , ...
%      'RearrangeableColumn','on', 'ColumnEditable', true, ...
%      'ToolTipString', [ char(a) NL 'You may edit the Signal from the object.' ...
%         NL 'The Contextual menu enables to plot a selection of the oject' ], ...    
    
    % attach contextual menu
    uicm = uicontextmenu; 
    S=a.Source; T=a.Title;
    if length(T) > 23, T=[ T(1:20) '...' ]; end
    if length(S) > 23, S=[ '...' S(end-20:end) ]; end
    uimenu(uicm, 'Separator','on', 'Label', [ 'Title: "' T '"' ]);
    uimenu(uicm, 'Label', [ 'Source: <' S '>' ]);
    uimenu(uicm, 'Label', [ 'Tag: '  a.Tag  ]);
    uimenu(uicm, 'Label', [ 'User: ' a.User ]);
    if ~isempty(a.Label)
      uimenu(uicm, 'Label', [ 'Label: ' a.Label ]);
    end
    if length(getaxis(a))
      uimenu(uicm, 'Separator','on', 'Label', '[Rank]         [Value] [Description]');
      for index=0:length(getaxis(a))
        [v, l] = getaxis(a, num2str(index));
        x      = getaxis(a, index);
        m      = get(a, 'Monitor');
        if index==0 & not(all(m==1 | m==0))
          uimenu(uicm, 'Label', sprintf('%6i %15s  %s [%g:%g] (per monitor)', index, v, l, min(x(:)), max(x(:))));
        else
          uimenu(uicm, 'Label', sprintf('%6i %15s  %s [%g:%g]', index, v, l, min(x(:)), max(x(:))));
        end
      end
    end
    uimenu(uicm, 'Separator','on','Label','Plot selection...', 'Callback', @iData_edit_display);
    uimenu(uicm, 'Separator','on','Label', 'About iData', 'Callback',[ 'msgbox(''' version(iData) ''')' ]);
    % attach contexual menu to the table
    try
      set(h,   'UIContextMenu', uicm); 
    end
  else
    % open the data file
    try
      edit(a.Source);
    end
  end
end

% ------------------------------------------------------------------------------
% Events
% ------------------------------------------------------------------------------

% function called when the Table selection is plotted
function iData_edit_display(hObject, eventdata, varargin)
  % hObject is the handle of the context menu item
  fig = gcbf;
  UserData  = get(fig, 'UserData');
  Signal    = get(UserData.handle, 'Data');
  Selection = UserData.Selection; % as { rows, columns }
  
  % get the Selection as Signal( rows, columns )
  if isempty(Selection) || isempty(Selection{1}) || isempty(Selection{2})
    Selection = Signal;
    rows    = [ 1 size(Signal,1) ];
    columns = [ 1 size(Signal,2) ];
  else
    rows    = Selection{1};
    columns = Selection{2};
    Selection = Signal(Selection{:});
  end
  clear Signal

  % plot it (no axes here ?)
  if isscalar(Selection) || ndims(Selection) > 2
    return;
  end
  nfig = figure;
  if isvector(Selection)
    plot(Selection);
    if numel(rows) == 1
      xlabel(sprintf('(%i,[%i:%i])', rows, min(columns),max(columns)));
    else
      xlabel(sprintf('([%i:%i],%i)', min(rows),max(rows), columns));
    end
  elseif ndims(Selection) == 2
    surf(Selection);
    ylabel(sprintf('[%i:%i]', min(rows),max(rows)));
    xlabel(sprintf('[%i:%i]', min(columns),max(columns)));
  end
  t = sprintf('Edit ([%i:%i,%i:%i]) %s', ...
        min(rows),max(rows),min(columns),max(columns), UserData.Name);
  set(nfig, 'Name', t);
  title(t)
  
end

function iData_edit_CellSelectionCallback(source, eventdata)
  % copy the selection into the UserData
  
  fig = gcbf; % current UItable figure
  
  rows    = unique(eventdata.Indices(:,1)); % rows and columns selected
  columns = unique(eventdata.Indices(:,2));
  
  UserData = get(fig, 'UserData');
  UserData.Selection = { rows,columns };
  set(fig, 'UserData',UserData);
end
