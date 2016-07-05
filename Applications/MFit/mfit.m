function varargout = mfit(varargin)
% MFIT M-file for mfit.fig
%      MFIT, by itself, creates a new MFIT or raises the existing
%      singleton*.
%
%      H = MFIT returns the handle to a new MFIT or the handle to
%      the existing singleton*.
%
%      MFIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MFIT.M with the given input arguments.
%
%      MFIT('Property','Value',...) creates a new MFIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mfit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mfit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mfit

% Last Modified by GUIDE v2.5 04-Jul-2016 14:58:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mfit_OpeningFcn, ...
                   'gui_OutputFcn',  @mfit_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before mfit is made visible.
function mfit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mfit (see VARARGIN)

% Choose default command line output for mfit
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mfit wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mfit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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



function Button_Data_File_Callback(hObject, eventdata, handles)
% hObject    handle to Button_Data_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Button_Data_File as text
%        str2double(get(hObject,'String')) returns contents of Button_Data_File as a double


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
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
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
