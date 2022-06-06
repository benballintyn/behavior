function varargout = sort_licks_gui(varargin)
% SORT_LICKS_GUI MATLAB code for sort_licks_gui.fig
%      SORT_LICKS_GUI, by itself, creates a new SORT_LICKS_GUI or raises the existing
%      singleton*.
%
%      H = SORT_LICKS_GUI returns the handle to a new SORT_LICKS_GUI or the handle to
%      the existing singleton*.
%
%      SORT_LICKS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SORT_LICKS_GUI.M with the given input arguments.
%
%      SORT_LICKS_GUI('Property','Value',...) creates a new SORT_LICKS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sort_licks_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sort_licks_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sort_licks_gui

% Last Modified by GUIDE v2.5 30-Dec-2018 19:18:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sort_licks_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @sort_licks_gui_OutputFcn, ...
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


% --- Executes just before sort_licks_gui is made visible.
function sort_licks_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sort_licks_gui (see VARARGIN)

% Choose default command line output for sort_licks_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sort_licks_gui wait for user response (see UIRESUME)
% uiwait(handles.sort_licks_gui);


% --- Outputs from this function are returned to the command line.
function varargout = sort_licks_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in keep_button.
function keep_button_Callback(hObject, eventdata, handles)
% hObject    handle to keep_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in reject_button.
function reject_button_Callback(hObject, eventdata, handles)
% hObject    handle to reject_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in select_data_menu.
function select_data_menu_Callback(hObject, eventdata, handles)
% hObject    handle to select_data_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns select_data_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from select_data_menu


% --- Executes during object creation, after setting all properties.
function select_data_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_data_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_data_button.
function load_data_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_data_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
