function hIcon=iView_private_mouse(instance, action, object)
% iView_private_icon_mouse Handle mouse/keyboard events in window/icons
% action: down, drag, up

if isempty(action), return; end

obj.type  =get(gco,'type');  
if strcmp(obj.type, 'uicontrol')
  obj.type =get(gco,'style');  
end

obj.tag   =get(gco,'tag');
press.type=get(gcf,'SelectionType'); % normal, extend(+shift), alternate(+ctrl), open(double click)

if isempty(obj.type) | isempty(obj.tag), return; end

% only acts on buttons (but not the 'center' one from the sliders)
if ~isempty(strmatch(obj.tag, {'Slider_v','Slider_h','Center'},'exact'))
  return;
end
if strcmp(action, 'down')
  obj
press
action
  switch obj.type
  case {'togglebutton','checkbox'}
    % action on an icon (open, drag'n drop)
    if strcmp(press.type, 'open')
      iview(gcf,'open',obj.tag);
    elseif strcmp(press.type, 'normal') | strcmp(press.type, 'alt')
      % handle movement of icon
      set(gcf,'Pointer','Fleur');                                  % Cursor to fleur
      %set(gcf,'WindowButtonMotionFcn',[ 'iview(' num2str(gcf) ',''mouse_drag'',''' obj.tag ''');' ]); % Drag when moved
      set(gcf,'WindowButtonUpFcn',[ 'iview(' num2str(gcf) ',''mouse_up'',''' obj.tag ''');' ]);  % end of movement
      setappdata(gcf,'dragplot',gco);
    end
  case 'figure'
    % action on the figure background (select rectangle)
    % sets WindowButtonMotionFcn and WindowButtonUpFcn
  end
elseif strcmp(action, 'drag')
  % do nothing yet. We just wait for end of movement. Could display here a message as 'Dragging object(s) so and so'
    obj
press
action
elseif strcmp(action, 'up')
  obj
press
action
  f = get(0,'pointerwindow');           % determine which figure window it's over (target). gcf=original
  h = getappdata(gcf,'dragplot');       % Retrieve the handle to the dragged object.
  o = gcf;
  
  if f > 0 & strcmp(get(f, 'Tag'),'iView_instance') % If it's a valid iView figure window,
    Data0  = getappdata(gcf, 'Data');     % data in original window
    Data1  = getappdata(f, 'Data');       % data in target window
    object = findobj(Data0, 'Tag', get(h,'Tag'));          % find iData which has same Tag as pushbutton
    p = get(f,'currentpoint');          % then determine the current position (target),
    t = get(gcf,'selectiontype');     % determine whether the CTRL key was pressed too.
    if iscell(object), object = object{1}; end
    index = iView_private_icon_index(gcf, p); % index of target position in the iData array
    index0= strmatch(object.Tag, get(Data0,'Tag'),'exact');
    if strcmp(t,'alt') & f~=gcf               % If CTRL was pressed (or right button),
      % make a copy of object (drag+ctrl)
      object = copyobj(object);
      iView_private_icon(f, 'load', object);  % add object at the end of all items (in target)
      % copy object to the target position
      if index > 1,              p1 = Data1(1:(index-1)); 
      else 
        if length(Data1) == 1, p1=Data1; 
        else p1= [] ; end
      end
      if index <= length(Data1) & length(Data1)>1, p2 = Data1(index:end);   else p2 = []; end
      Data1 = [ p1 object p2];
      setappdata(f, 'Data', Data1);
      % re-order icons
      iView_private_icon(f, 'resize', []);
    else
      if f==gcf 
        % move in same figure
        % index of original object
        if length(Data0)>1 & index ~= index0
          index0= strmatch(object.Tag, get(Data0,'Tag'),'exact');
          Data0(index0) = [];  % remove original position
          if index > 1,              p1 = Data0(1:(index-1)); 
          else 
            if length(Data0) == 1, p1=Data0; 
            else p1= [] ; end
          end
          if index <= length(Data0) & length(Data0)>1, p2 = Data0(index:end);   else p2 = []; end
          Data0 = [ p1 object p2];
          setappdata(f, 'Data', Data0);
          % re-order icons
          iView_private_icon(f, 'resize', []);
        end
      else
        % get object, copy it to target window, remove original
        % make a copy of object (drag+ctrl)
        iView_private_icon(f, action, object);  % add object at the end of all items (in target)
        % copy object to the target position
        if index > 1,              p1 = Data1(1:(index-1)); 
        else 
          if length(Data1) == 1, p1=Data1; 
          else p1= [] ; end
        end
        if index <= length(Data1) & length(Data1)>1, p2 = Data1(index:end);   else p2 = []; end
        Data1 = [ p1 object p2];
        setappdata(f, 'Data', Data1);
        % re-order icons
        iView_private_icon(f, 'resize', []);  % target window is updated
        if length(Data0) > 1,
          Data0(index0) = [];  % remove original position
        else
          Data0=[];
        end
        setappdata(gcf, 'Data', Data0);
        iView_private_icon(f, 'resize', []);  % target window is updated
        iView_private_icon(f, 'check', []);   % target window is updated
      end
    end
  end
  
  set(gcf,'WindowButtonMotionFcn','');
  set(gcf,'Pointer','arrow');
  set(gcf,'WindowButtonUpFcn','');
end


function index = iView_private_icon_index(instance, position)
  config = iView_private_config(instance, 'load');

  % size of Icon is [ 2*config.IconSize config.IconSize ]

  % number of elements per row
  instancePosition = get(instance, 'Position');
  instanceWidth    = instancePosition(3);
  nbElementsPerRow = max(1,floor(instanceWidth/(2.25*config.IconSize)));
  
  % determine col/row indexes from position of pointer
  iColumn = position(1) - config.IconSize/4;
  iColumn = max(1,floor(iColumn / (2.25*config.IconSize)));
  iRow    = position(2) - config.IconSize/4;
  iRow    = max(1,floor(iRow/ (1.25*config.IconSize)));
  
  index   = ((iRow-1)*nbElementsPerRow+(iColumn-1))+1;
  index   = min(length(getappdata(instance, 'Data')), index);


