function varargout = initialize(varargin)
% INITIALIZE MATLAB code for initialize.fig
%      INITIALIZE, by itself, creates a new INITIALIZE or raises the existing
%      singleton*.
%
%      H = INITIALIZE returns the handle to a new INITIALIZE or the handle to
%      the existing singleton*.
%
%      INITIALIZE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INITIALIZE.M with the given input arguments.
%
%      INITIALIZE('Property','Value',...) creates a new INITIALIZE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before initialize_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to initialize_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help initialize

% Last Modified by GUIDE v2.5 19-Oct-2011 12:59:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @initialize_OpeningFcn, ...
                   'gui_OutputFcn',  @initialize_OutputFcn, ...
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

% --- Executes just before initialize is made visible.
function initialize_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to initialize (see VARARGIN)

% Choose default command line output for initialize
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using initialize.
if strcmp(get(hObject,'Visible'),'off')
axes(handles.axes1);
cla;
size  = str2num(get(handles.edit2, 'String')) ;
positions  = explode(get(handles.edit1, 'String'),';') ;
posarray = [] ; 
for i=1:1:numel(positions)
    posarray(i) = str2num(strtrim(positions{i})) ;
end ;
hist(posarray,size) ;
axis([0 size 0 1]) ;


u_min = str2num(get(handles.edit3, 'String')) ;
u_max= str2num(get(handles.edit10, 'String')) ;
velocity_u = [] ;
if(get(handles.radiobutton2, 'Value')==1)
   
    velocity_u = u_max-(0:(u_max-u_min)/(size):u_max-u_min) ;
    if(u_max==u_min)
        velocity_u = ones(1,size+1)*u_max ;
    end ;
elseif(get(handles.radiobutton1, 'Value')==1)
    velocity_u = (u_min:(u_max-u_min)/(size):u_max) ;
    if(u_max==u_min)
        velocity_u = ones(1,size+1)*u_max ;
    end 
elseif(get(handles.radiobutton3, 'Value')==1)
    mid = floor((size+1)/2) ;
    a = (u_max-u_min)/(mid^2) ;
    velocity_u = u_min+a*(-mid:1:-mid+size).^2 ;
end
axes(handles.axes2);
cla;
plot(velocity_u) ;

end

% UIWAIT makes initialize wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = initialize_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;
size  = str2num(get(handles.edit2, 'String')) ;
positions  = explode(get(handles.edit1, 'String'),';') ;
posarray = [] ; 
for i=1:1:numel(positions)
    posarray(i) = str2num(strtrim(positions{i})) ;
end ;
u = str2num(get(handles.edit3, 'String')) ;
dx = str2num(get(handles.edit4, 'String')) ;
dt = str2num(get(handles.edit5, 'String')) ;
sc = str2num(get(handles.edit6, 'String')) ;
it = str2num(get(handles.edit7, 'String')) ;
N = str2num(get(handles.edit9, 'String')) ;
fid = fopen(get(handles.edit8, 'String'), 'w');
fprintf(fid, 'M %d\n', size);
fprintf(fid, 'N %d\n', N);
fprintf(fid, 'dt %f\n', dt);
fprintf(fid, 'dx %f\n', dx);
fprintf(fid, 'eps %f\n', 0.000001);
fprintf(fid, 'uv %f\n', u);
fprintf(fid, 'sc %f\n', sc);
fprintf(fid, 'iterations %d\n', it);
initial = zeros(1,size) ; 
for i = 1 : numel(posarray) ;
    initial(posarray(i)) =  1 ;
end ;
fprintf(fid, 'initial_x');
for i = 1 : size ; 
    fprintf(fid, ' %d', initial(i));
end ;
fclose(fid);

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
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
size  = str2num(get(handles.edit2, 'String')) ;
positions  = explode(get(handles.edit1, 'String'),';') ;
posarray = [] ; 
for i=1:1:numel(positions)
    posarray(i) = str2num(strtrim(positions{i})) ;
end ;
hist(posarray,size) ;
axis([0 size 0 1]) ;


u_min = str2num(get(handles.edit3, 'String')) ;
u_max= str2num(get(handles.edit10, 'String')) ;
velocity_u = [] ;
if(get(handles.radiobutton2, 'Value')==1)
   
    velocity_u = u_max-(0:(u_max-u_min)/(size):u_max-u_min) ;
    if(u_max==u_min)
        velocity_u = ones(1,size+1)*u_max ;
    end ;
elseif(get(handles.radiobutton1, 'Value')==1)
    velocity_u = (u_min:(u_max-u_min)/(size):u_max) ;
    if(u_max==u_min)
        velocity_u = ones(1,size+1)*u_max ;
    end 
elseif(get(handles.radiobutton3, 'Value')==1)
    mid = floor((size+1)/2) ;
    a = (u_max-u_min)/(mid^2) ;
    velocity_u = u_min+a*(-mid:1:-mid+size).^2 ;
end
title('Initial situation psi_0') ; 
xlabel('Gridbox i') ;
ylabel('psi_0') ;
axes(handles.axes2);
cla;
plot(velocity_u) ;
title('Velocity vector u') ; 
xlabel('Gridbox border i+1/2') ;
ylabel('u_i+1/2') ;
dx = str2num(get(handles.edit4, 'String')) ;
dt = str2num(get(handles.edit5, 'String')) ;
sc = str2num(get(handles.edit6, 'String')) ;
it = str2num(get(handles.edit7, 'String')) ;
if(max((abs(velocity_u)*dt)/dx)>1)
    set(handles.text8, 'String','UNSTABLE') ;
    set(handles.text8, 'ForegroundColor','red') ;

else
    set(handles.text8, 'String','STABLE') ;
    set(handles.text8, 'ForegroundColor','green') ;
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton4.
function pushbutton4_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
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
size  = str2num(get(handles.edit2, 'String')) ;
positions  = explode(get(handles.edit1, 'String'),';') ;
posarray = [] ; 
for i=1:1:numel(positions)
    posarray(i) = str2num(strtrim(positions{i})) ;
end ;
u_min = str2num(get(handles.edit3, 'String')) ;
u_max= str2num(get(handles.edit10, 'String')) ;
dx = str2num(get(handles.edit4, 'String')) ;
dt = str2num(get(handles.edit5, 'String')) ;
sc = str2num(get(handles.edit6, 'String')) ;
it = str2num(get(handles.edit7, 'String')) ;
N = str2num(get(handles.edit9, 'String')) ;
velocity_u = [] ;
if(get(handles.radiobutton2, 'Value')==1)

    velocity_u = u_max-(0:(u_max-u_min)/(size):u_max-u_min) ;
    if(u_max==u_min)
        velocity_u = ones(1,size+1)*u_max ;
    end ;
elseif(get(handles.radiobutton1, 'Value')==1)
    velocity_u = (u_min:(u_max-u_min)/(size):u_max) ;
    if(u_max==u_min)
        velocity_u = ones(1,size+1)*u_max ;
    end 
elseif(get(handles.radiobutton3, 'Value')==1)
    mid = floor((size+1)/2) ;
    a = (u_max-u_min)/(mid^2) ;
    velocity_u = u_min+a*(-mid:1:-mid+size).^2 ;
end

% Write first file
fid = fopen(get(handles.edit8, 'String'), 'w');
fprintf(fid, 'M %d\n', size);
fprintf(fid, 'N %d\n', N);
fprintf(fid, 'dt %f\n', dt);
fprintf(fid, 'dx %f\n', dx);
fprintf(fid, 'eps %f\n', 0.000001);
fprintf(fid, 'sc %f\n', sc);
fprintf(fid, 'iterations %d\n', it);
initial = zeros(1,size) ; 
for i = 1 : numel(posarray) ;
    initial(posarray(i)) =  1 ;
end ;
fprintf(fid, 'initial_x');
for i = 1 : size ; 
    fprintf(fid, ' %d', initial(i));
end ;
fprintf(fid,'\n') ;

disp(velocity_u) ;

fprintf(fid, 'velocity_u');
for i = 1 : size+1 ; 
    fprintf(fid, ' %f',velocity_u(i));
end ;


fclose(fid);




% Write second file
fid = fopen('noit.txt', 'w');
fprintf(fid, 'M %d\n', size);
fprintf(fid, 'N %d\n', N);
fprintf(fid, 'dt %f\n', dt);
fprintf(fid, 'dx %f\n', dx);
fprintf(fid, 'eps %f\n', 0.000001);
fprintf(fid, 'sc %f\n', sc);
fprintf(fid, 'iterations %d\n', 0);
initial = zeros(1,size) ; 
for i = 1 : numel(posarray) ;
    initial(posarray(i)) =  1 ;
end ;
fprintf(fid, 'initial_x');
for i = 1 : size ; 
    fprintf(fid, ' %d', initial(i));
end ;
fprintf(fid,'\n') ;


fprintf(fid, 'velocity_u');
for i = 1 : size+1 ; 
    fprintf(fid, ' %f',velocity_u(i));
end ;


fclose(fid);

stat = unix(['../smolarki.out ' get(handles.edit8, 'String')]) ;
stat = unix('../smolarki.out noit.txt') ;

disp(stat) ;
wave = rand(100,240) ;
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

figure  ;
for i = 1 : N ; 
    clf ;
    plot(wave0(i,:)) ;  
    hold on ;
    plot(mywave(i,:),'r') ;
    axis([0 size 0 1]) ;
    legend('No antidiffusion', ['Antidiffusion with ' num2str(it) ' iterations']) ;
    drawnow ; 
end ;



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4
