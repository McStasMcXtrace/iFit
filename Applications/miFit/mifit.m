function out = mifit(varargin)

    out = [];

    % look if the main window is opened
    fig = mifit_fig;
    if isempty(fig) || ~ishandle(fig)
        fig = feval([ mfilename '_OpeningFcn' ]);
        out = fig;
    end

    if nargin == 1 && isa(varargin{1}, 'iFunc')
        % import model: add it to the Models menu, if not already there
    elseif nargin == 1 && isa(varargin{1}, 'iData')
        % import iData object(s) into the List
        disp('Import into List')
    end

    if ~isempty(varargin) && ischar(varargin{1})
        % callback with varargin{1} == 'action'
        if strcmpi(varargin{1},'identify'), return; end
        try
          feval(varargin{:});
        catch
          disp([ mfilename ': Unknown action "' varargin{1} '"' ])
        end
    end
    
% -------------------------------------------------------------------------
% Main window ID and creation
% -------------------------------------------------------------------------

function f=mifit_fig(tag)
% search for a given Tag in Application or main Figure if ommitted.
  persistent fig handles
  
  if ~ishandle(fig), fig=[]; end
  if isempty(fig)
    fig = findall(0, 'Tag','MiFit');
    if length(fig) > 1, delete(fig(2:end)); end % unique instance
    handles = [];
  end

  if nargin == 0
    f=fig;
  else
    if ~isfield(handles, tag) || ~ishandle(handles.(tag))
      handles.(tag) = findobj(fig, 'Tag', tag);
    end
    f = handles.(tag);
  end

% --- Executes just before mifit is made visible.
function fig = mifit_OpeningFcn
% This function creates the main window and returns its ID

fig = mifit_fig;
if isempty(fig) || ~ishandle(fig)
    % create the main figure
    disp([ mfilename ': Welcome to MiFit !' ])
    disp(datestr(now))
    fig = openfig(mfilename);
    
    % test if any of the Models/Optimizers menu is empty
    % if not, return (no need to build them)
    hmodels = mifit_fig('Menu_Models');
    hoptim  = mifit_fig('Menu_Optimizers');
    if isempty(hmodels) || isempty(hoptim) ...
            || ~isempty([ get(hmodels,'Children') ; get(hoptim,'Children') ]), return; end
    
    % Display welcome dialog during menu build
    h = msgbox('Welcome to MiFit !','MiFit: Starting','help');
    
    % get the list of Models and Optimizers
    [optimizers,functions] = fits(iFunc);
    
    % fill Models menu
    if any(~isempty(functions)) && all(isa(functions, 'iFunc'))
        disp([ mfilename ': Initializing Models...' ]);
        separator = 1;
        % first add the 'Add new Model': from File or ifitmakefunc
        uimenu(hmodels, 'Label','Add new Model...');
        for f=functions
            % each Model is an iFunc object. These should be stored in the
            % Models items 'UserData'
            if ~isempty(f) && ~isempty(f.Name)
                if separator
                  uimenu(hmodels, 'Label', f.Name, 'UserData', f, 'Separator','on');
                  separator = 0;
                else
                  uimenu(hmodels, 'Label', f.Name, 'UserData', f);
                end
            end
        end
    end
    
    % fill Optimizers menu
    if ~isempty(optimizers) && iscell(optimizers)
        disp([ mfilename ': Initializing Optimizers...' ]);
        for f=optimizers
            % each optimizer is given with its function name. We request
            % 'defaults' and display its name
            
            o=feval(f{1},'defaults');
            if isfield(o, 'algorithm') && ~isempty(o.algorithm)
                algorithm = o.algorithm;
            else
                algorithm = f{1};
            end
            if ~isempty(algorithm)
                uimenu(hoptim, 'Label', algorithm, 'UserData', f{1});
            end
        end
    end
    
    % close welcome image
    delete(h);
end

% -------------------------------------------------------------------------
% Callbacks
% -------------------------------------------------------------------------

% --- Executes on button press in Button_Fit.
function Button_Fit_Callback(hObject, eventdata, handles)
% hObject    handle to Button_Fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Button_Plot.
function Button_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to Button_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function Button_Data_File_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Button_Data_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Menu_File_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Menu_Load_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Menu_Models_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_Models (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Menu_Optimizers_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_Optimizers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Menu_Tools_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_Tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Tools_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to Tools_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Tools_Fit_Callback(hObject, eventdata, handles)
% hObject    handle to Tools_Fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function File_New_Callback(hObject, eventdata, handles)
% hObject    handle to File_New (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function File_Open_Callback(hObject, eventdata, handles)
% hObject    handle to File_Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  List_Add_Callback;


% --------------------------------------------------------------------
function File_Save_Callback(hObject, eventdata, handles)
% hObject    handle to File_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function File_Preferences_Callback(hObject, eventdata, handles)
% hObject    handle to File_Preferences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function File_Exit_Callback(hObject, eventdata, handles)
% hObject    handle to File_Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  disp([ mfilename ': Exiting MiFit, saving Preferences and Data files...' ])
  disp(datestr(now))


% --- Executes on selection change in List_Data_Files.
function List_Data_Files_Callback(hObject, eventdata, handles)
% hObject    handle to List_Data_Files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns List_Data_Files contents as cell array
%        contents{get(hObject,'Value')} returns selected item from List_Data_Files


% --- Executes during object creation, after setting all properties.
function List_Data_Files_CreateFcn(hObject, eventdata, handles)
% hObject    handle to List_Data_Files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Button_Open.
function Button_Open_Callback(hObject, eventdata, handles)
% hObject    handle to Button_Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  List_Add_Callback;


% --- Executes on button press in Button_Data_View.
function Button_Data_View_Callback(hObject, eventdata, handles)
% hObject    handle to Button_Data_View (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function List_Add_Callback(hObject, eventdata, handles)
% hObject    handle to List_Add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  disp([ mfilename ': Add Data files to the MiFit List' ]);
  toadd = iData('');
  mifit(toadd);

% --------------------------------------------------------------------
function List_Remove_Callback(hObject, eventdata, handles)
% hObject    handle to List_Remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function List_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to List_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function List_Fit_Callback(hObject, eventdata, handles)
% hObject    handle to List_Fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Tools_View_Callback(hObject, eventdata, handles)
% hObject    handle to Tools_View (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Tools_Help_Callback(hObject, eventdata, handles)
% hObject    handle to Tools_Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Tools_About_Callback(hObject, eventdata, handles)
% hObject    handle to Tools_About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    h = msgbox('MiFit is a user interface to the iFit software','MiFit: About','help');

% --------------------------------------------------------------------
function File_Close_Callback(hObject, eventdata, handles)
% hObject    handle to File_Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Menu_List_Data_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_List_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function List_Remove_All_Callback(hObject, eventdata, handles)
% hObject    handle to List_Remove_All (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function List_Select_All_Callback(hObject, eventdata, handles)
% hObject    handle to List_Select_All (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
