function hIcon=iView_private_mouse(instance, action, object)
% iView_private_icon_mouse Handle mouse/keyboard events in window/icons
% action: down, drag, up

if isempty(action), return; end

config = iView_private_config(instance, 'load');

obj.handle= gco;
obj.type  =get(gco,'type');  
if strcmp(obj.type, 'uicontrol')
  obj.type =get(gco,'style');  
end

obj.Tag   =get(gco,'tag');
% event: normal, extend(+shift), alternate(+ctrl), open(double click)
obj.event=get(gcf,'SelectionType'); 

if isempty(obj.type) | isempty(obj.Tag), return; end

% only acts on buttons (but not the 'center' one from the sliders)
if ~isempty(strmatch(obj.Tag, {'Slider_v','Slider_h','Center'},'exact'))
  return;
end

switch action
case 'down'
  set(gcf,'WindowButtonMotionFcn','');
  set(gcf,'WindowButtonUpFcn','');
  set(gcf,'Pointer','arrow');
  set(obj.handle,'Selected','off');
  switch obj.type
  case {'togglebutton','checkbox','frame','pushbutton','text'}
    % action on an icon (open, drag'n drop)
    set(obj.handle, 'Enable', 'inactive');
    index=iView_private_icon_index(gcf, get(gcf, 'CurrentPoint'));
    if strcmp(obj.event, 'open') | strcmp(obj.event, 'normal')
      % an open is issued after a down+up
      % counter balance initial Down
      v = get(obj.handle, 'Value');
      if v==0, v=1; else v=0; end
      set(obj.handle, 'Value', v);
      % open ?
      if strcmp(obj.event, 'open')
        iview(gcf,'data_open','selection');
        % open the other selected objects
      end
    end
    if strcmp(obj.event, 'normal') | strcmp(obj.event, 'alt')
      % handle movement of icon
      set(gcf,'WindowButtonUpFcn',[ 'iview(' num2str(gcbf) ',''mouse_up'',''' obj.Tag ''');' ]);  % end of movement
      set(gcf,'WindowButtonMotionFcn',[ 'iview(' num2str(gcbf) ',''mouse_drag'',''' obj.Tag ''');' ]);
      obj.pointer = get(gcf, 'CurrentPoint');
      obj.dragged = 0;
      setappdata(gcf, 'iView_motion_obj',obj);
    end
  case 'figure'
    if strcmp(obj.event, 'normal')
      % action on the figure background (select rectangle)
      finalbox = rbbox;
      % find objects within area and toggle their Value
      objects = get(gcf,'Children');
      for index=1:length(objects)
        this=objects(index);
        if strcmp(get(this,'type'),'uicontrol')
          if ~isempty(strmatch(get(this,'Style'), {'togglebutton','checkbox','pushbutton'}))
            pos = get(this,'Position');
            if all(pos(1:2) > finalbox(1:2)) & all(pos(1:2)+pos(3:4) < finalbox(1:2)+finalbox(3:4))
              v = get(this, 'Value');
              if v==0, v=1; else v=0; end
              set(this, 'Value', v);
            end
          end
        end
      end
    end
  end
  
case 'drag'
  % We just wait for end of movement. Could display here a message as 'Dragging object(s) so and so'
  if ~isempty(strmatch(obj.type, {'togglebutton','checkbox','frame','pushbutton','text'}))
    % check if we have moved enought from intial position
    obj0=getappdata(gcf,'iView_motion_obj');
    if any(abs(get(gcf,'CurrentPoint') - obj0.pointer) >= 5)
      % unactivate UIContextMenu for movement
      cmenu = get(obj.handle, 'UIContextMenu');
      set(cmenu, 'Visible', 'off');
      % set up pointer for movement
      set(gcf,'WindowButtonMotionFcn',''); % once we know we move, unactivate event
      if strcmp(get(gcf,'Pointer'),'arrow')
        if strcmp(obj.event, 'normal'), set(gcf,'Pointer','Fleur'); % Cursor to fleur for moving
        else set(gcf,'Pointer','crosshair'); end  % Cross pointer when copying
      end
      obj0.dragged=1;
      % highlight moved item
      set(obj.handle,'Selected','on', 'Value', 1);
    end
    setappdata(gcf, 'iView_motion_obj',obj0);
  end
  
case 'up'
  fig_f = get(0,'pointerwindow');           % determine which figure window it's over (target).
  fig_i = gcf;                              % gcf=original figure
  % obj contains object which is dragged
  set(fig_i,'WindowButtonMotionFcn','');
  set(fig_i,'WindowButtonUpFcn','');
  set(fig_i,'Pointer','arrow');
  set(obj.handle,'Selected','off', 'Enable', 'inactive');
  
  obj0=getappdata(gcf,'iView_motion_obj'); % information on original event
  
  % check if object was dragged
  if obj0.dragged
    % unactivate UIContextMenu for movement
    cmenu = get(obj.handle, 'UIContextMenu');
    set(cmenu, 'Visible', 'off');
  end

  if obj0.dragged && strcmp(get(fig_f,'Tag'), 'iView_instance') ...
    && ~isempty(strmatch(obj.type, {'togglebutton','checkbox','frame','pushbutton','text'}))
    % deposit object in Figure fig_f
    % position of target
    finalpos = get(fig_f, 'CurrentPoint');
    targetIndex=iView_private_icon_index(fig_f, finalpos);
    % check if we drop on background or on an object
    obj0.onbackground=1; hIcon=[];
    
    % get iData(targetIndex).Tag, Position (in target window)
    Data_i = getappdata(fig_i, 'Data');
    Data_f = getappdata(fig_f, 'Data');
    if targetIndex>0 & targetIndex<=length(Data_f)
      object = ind2sub(Data_f, targetIndex);
      % object is the Data item below cursor
      if isa(object, 'iData')
        hIcon = findobj(fig_f, 'Tag', object.Tag);
        pos   = get(hIcon, 'Position');
        % cursor on indexed object in target ?
        if all(pos(1:2) < finalpos(1:2)) & all(pos(1:2)+pos(3:4) > finalpos(1:2))
          obj0.onbackground=0;
        end
      end
    end

    % get initial selection (all Value=1 objects).
    [selection,selectedIndex,selected]= iView_private_selection(fig_i); % selected iData objects (initial)
    
    if isempty(selection) return; end
    
    % in the case of a move, we must remove selected items
    if obj0.onbackground && strcmp(obj.event, 'normal')
      iView_private_icon(fig_i, 'delete', selection);
      Data_i = getappdata(fig_i, 'Data');
      if fig_i == fig_f
        Data_f = Data_i;
      end
    end
    target_before_index = ind2sub(Data_f, 1:(targetIndex-1));          % iData before target position
    target_after_index  = ind2sub(Data_f, targetIndex:length(Data_f)); % iData after target position
    
    if ~obj0.onbackground
      % to object: merge/combine (normal) or custom (alt) (cat, overlay plot, minus, plus, ...)
      % get target iData from targetIndex
      % do operation
      % create resulting iData and position next to target
    else
      % to background: move or copy ?
      if strcmp(obj.event, 'alt')
        % copy
        selection = copyobj(selection);
      end
      % move/copy
      Data_f = [ target_before_index selection target_after_index ];
      setappdata(fig_f, 'Data', Data_f);
      iView_private_documents(fig_f);
      [hIcon,config,Data_f]=iView_private_icon(fig_f, 'check_load', Data_f);
    end
  end % if obj0.dragged
  
end % switch action

% =============================================================================
%                    INLINE functions
% =============================================================================
% iView_private_icon_index
% iView_private_icon_reorder
% iView_private_icon_selected

function index = iView_private_icon_index(instance, position)
% iView_private_icon_index Determine index of object below cursor position
% returned index is within 0:length(array)+1

  config = iView_private_config(instance, 'load');

  % size of Icon is [ config.IconWidth*config.IconSize config.IconSize ]

  % number of elements per row
  instancePosition = get(instance, 'Position');
  instanceWidth    = instancePosition(3);
  % reverse Height
  position(2) = instancePosition(4) - position(2);
  nbElementsPerRow = max(1,floor(instanceWidth/((config.IconWidth+.25)*config.IconSize)));
  
  % determine col/row indexes from position of pointer
  iColumn = position(1) - config.IconSize/4;
  iColumn = floor(iColumn / ((config.IconWidth+.25)*config.IconSize));
  iRow    = position(2) + config.IconSize/4;
  iRow    = floor(iRow/ (1.25*config.IconSize));
  
  index   = iRow*nbElementsPerRow+iColumn+1;
  index   = min(length(getappdata(instance, 'Data'))+1, index);


