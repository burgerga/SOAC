function varargout = initialized2d(varargin)
% INITIALIZED2D MATLAB code for initialized2d.fig
%      INITIALIZED2D, by itself, creates a new INITIALIZED2D or raises the existing
%      singleton*.
%
%      H = INITIALIZED2D returns the handle to a new INITIALIZED2D or the handle to
%      the existing singleton*.
%
%      INITIALIZED2D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INITIALIZED2D.M with the given input arguments.
%
%      INITIALIZED2D('Property','Value',...) creates a new INITIALIZED2D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before initialized2d_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to initialized2d_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help initialized2d

% Last Modified by GUIDE v2.5 27-Oct-2011 14:46:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @initialized2d_OpeningFcn, ...
                   'gui_OutputFcn',  @initialized2d_OutputFcn, ...
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


% --- Executes just before initialized2d is made visible.
function initialized2d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to initialized2d (see VARARGIN)

% Choose default command line output for initialized2d
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes initialized2d wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = initialized2d_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cone = get(handles.checkbox1, 'Value') 


size  = str2num(get(handles.edit16, 'String')) ;
position = size*(size/2)+(size/2) ;
csize  = str2num(get(handles.edit13, 'String')) ;
u = str2num(get(handles.edit14, 'String')) ;
v = str2num(get(handles.edit15, 'String')) ;
dx = str2num(get(handles.edit17, 'String')) ;
dy = str2num(get(handles.edit18, 'String')) ;
dt = str2num(get(handles.edit19, 'String')) ;
sc = str2num(get(handles.edit22, 'String')) ;
it = str2num(get(handles.edit20, 'String')) ;
N = str2num(get(handles.edit21, 'String')) ;

% Write first file
fid = fopen('input.txt', 'w');
fprintf(fid, 'M %d\n', size);
fprintf(fid, 'N %d\n', N);
fprintf(fid, 'dt %f\n', dt);
fprintf(fid, 'dx %f\n', dx);
fprintf(fid, 'dy %f\n', dy);
fprintf(fid, 'eps %f\n', 0.000001);
fprintf(fid, 'sc %f\n', sc);
fprintf(fid, 'cone %d\n', cone);
fprintf(fid, 'u %f\n', u);
fprintf(fid, 'v %f\n', v);
fprintf(fid, 'iterations %d\n', it);
fprintf(fid, 'initial_pos %d\n', position) ;
fprintf(fid, 'cloudsize %d\n', csize) ;

fclose(fid);




% Write second file
fid = fopen('noit.txt', 'w');
fprintf(fid, 'M %d\n', size);
fprintf(fid, 'N %d\n', N);
fprintf(fid, 'dt %f\n', dt);
fprintf(fid, 'dx %f\n', dx);
fprintf(fid, 'dy %f\n', dy);
fprintf(fid, 'eps %f\n', 0.000001);
fprintf(fid, 'sc %f\n', sc);
fprintf(fid, 'u %f\n', u);
fprintf(fid, 'cone %d\n', cone);
fprintf(fid, 'v %f\n', v);
fprintf(fid, 'iterations %d\n', 0);
fprintf(fid, 'initial_pos %d\n', position) ;
fprintf(fid, 'cloudsize %d\n', csize) ;

fclose(fid);


stat = unix('../smolarki2d.out input.txt') ;
stat = unix('../smolarki2d.out noit.txt') ;
load('wave0.dat') ;
mywave = [] ;
switch it
   case 1
      load('wave1.dat') ;
      mywave = wave1 ;
   case 2
      load('wave2.dat') ;
      mywave = wave2 ;
   case 3
      load('wave3.dat') ;
      mywave = wave3 ;
   otherwise
      disp('Three iterations max.')
end
maxx = max(max(mywave)) ;
for i = 1 : N-1 ; 
axes(handles.axes3);
cla;
contour(reshape(wave0(i,:),size,size)) ; caxis([0 maxx]) ; colorbar;
axes(handles.axes4);
cla;
contour(reshape(mywave(i,:),size,size)) ; caxis([0 maxx]) ; colorbar;
end ;
function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
u = str2num(get(handles.edit14, 'String')) ;
v = str2num(get(handles.edit15, 'String')) ;
dx = str2num(get(handles.edit17, 'String')) ;
dy = str2num(get(handles.edit18, 'String')) ;
dt = str2num(get(handles.edit19, 'String')) ;
sc = str2num(get(handles.edit22, 'String')) ;
it = str2num(get(handles.edit20, 'String')) ;
N = str2num(get(handles.edit21, 'String')) ;

if(sqrt((u^2*dt^2)/(dx^2)+(v^2*dt^2)/(dy^2))>(2^(-1/2)))
    set(handles.text22, 'String','UNSTABLE') ;
    set(handles.text22, 'ForegroundColor','red') ;

else
    set(handles.text22, 'String','STABLE') ;
    set(handles.text22, 'ForegroundColor','green') ;
end



% --- Executes during object creation, after setting all properties.
function text22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
it = str2num(get(handles.edit20, 'String')) ;
N = str2num(get(handles.edit21, 'String')) ;
size  = str2num(get(handles.edit16, 'String')) ;

load('wave0.dat') ;
mywave = [] ;
switch it
   case 1
      load('wave1.dat') ;
      mywave = wave1 ;
   case 2
      load('wave2.dat') ;
      mywave = wave2 ;
   case 3
      load('wave3.dat') ;
      mywave = wave3 ;
   otherwise
      disp('Three iterations max.')
end

for i = 1 : N-1 ; 
axes(handles.axes3);
cla;
contour(reshape(wave0(i,:),size,size),50) ; colorbar;
axes(handles.axes4);
cla;
contour(reshape(mywave(i,:),size,size),50) ; colorbar;
end ;
