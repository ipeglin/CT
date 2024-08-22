function varargout = ART_GUI(varargin)
% ART_GUI MATLAB code for ART_GUI.fig
%      ART_GUI, by itself, creates a new ART_GUI or raises the existing
%      singleton*.
%
%      H = ART_GUI returns the handle to a new ART_GUI or the handle to
%      the existing singleton*.
%
%      ART_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ART_GUI.M with the given input arguments.
%
%      ART_GUI('Property','Value',...) creates a new ART_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ART_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ART_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ART_GUI

% Last Modified by GUIDE v2.5 02-Dec-2017 12:23:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ART_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ART_GUI_OutputFcn, ...
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


% --- Executes just before ART_GUI is made visible.
function ART_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ART_GUI (see VARARGIN)

% Choose default command line output for ART_GUI
handles.output = hObject;
set(handles.size,'Enable','off');
set(handles.iteration,'Enable','off');
set(handles.lambda,'Enable','off');
set(handles.alpha,'Enable','off');
set(handles.snr,'Enable','off');
set(handles.reconstruct,'Enable','off');
% Update handles structure

guidata(hObject, handles);

% UIWAIT makes ART_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ART_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in import.
function import_Callback(hObject, eventdata, handles)
% hObject    handle to import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname] = uigetfile('*.mat','Select projection data file');

if isequal(filename,0)||isequal(pathname,0)
    return
else
    projection_data=importdata(filename);
end

axes(handles.sinogram);
imagesc(projection_data);
title('Sinogram');
colorbar;
xlabel('projection angle');
ylabel('ray position')
handles.projection_data=projection_data;
set(handles.size,'Enable','on');
guidata(hObject,handles);





function size_Callback(hObject, eventdata, handles)
% hObject    handle to size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val2=get(handles.size,'Value');

switch (val2)
    case 1
        image_size=128;
        ray_number=183;
    case 2
        image_size=50;
        ray_number=73;
end

projection_data=Func_sinogram(handles.projection_data,image_size);


set(handles.iteration,'Enable','on');



handles.image_size=image_size;
handles.projection_data=projection_data;
handles.ray_number=ray_number;
guidata(hObject,handles);


% Hints: get(hObject,'String') returns contents of size as text
%        str2double(get(hObject,'String')) returns contents of size as a double


% --- Executes during object creation, after setting all properties.
function size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iteration_Callback(hObject, eventdata, handles)
% hObject    handle to iteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

iteration_number=str2num(get(hObject,'String'));
set(handles.lambda,'Enable','on');

handles.iteration_number=iteration_number;
guidata(hObject,handles);

% Hints: get(hObject,'String') returns contents of iteration as text
%        str2double(get(hObject,'String')) returns contents of iteration as a double


% --- Executes during object creation, after setting all properties.
function iteration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lambda_Callback(hObject, eventdata, handles)
% hObject    handle to lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
relaxation_factor=str2double(get(hObject,'String'));
set(handles.alpha,'Enable','on');

% Hints: get(hObject,'String') returns contents of lambda as text
%        str2double(get(hObject,'String')) returns contents of lambda as a double
handles.relaxation_factor=relaxation_factor;
guidata(hObject,handles);


function snr_Callback(hObject, eventdata, handles)
% hObject    handle to snr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
noise_level=str2double(get(hObject,'String'));

set(handles.reconstruct,'Enable','on');
handles.noise_level=noise_level;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in alpha.
function alpha_Callback(hObject, eventdata, handles)
set(handles.snr,'Enable','on');
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val1=get(handles.alpha,'Value');

switch (val1)
    case 1
%         w=w_01;
        w_name='discrete model';
    case 2
%         w=w_length;
        w_name='line based model';
    case 3
%         w=w_area;
        w_name='area based model';
end



handles.val1=val1;
handles.w_name=w_name;
% handles.w=w;
guidata(hObject,handles);

% Hints: contents = cellstr(get(hObject,'String')) returns alpha contents as cell array
%        contents{get(hObject,'Value')} returns selected item from alpha


% --- Executes during object creation, after setting all properties.
function alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: size controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in reconstruct.
function reconstruct_Callback(hObject, eventdata, handles)


% hObject    handle to reconstruct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = waitbar(0,'loading parameters','Name','',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)

if handles.image_size==128
    
load('WeightingMatrix.mat');
else if handles.image_size==50
    load ('WeightingMatrix50.mat');
    end
end

switch (handles.val1)
    case 1
        w=w_01;
        
    case 2
        w=w_length;
        
    case 3
        w=w_area;
        
end






projection_number=size(handles.projection_data,2);

[weighting_factor]=build_w(w,projection_number,handles.image_size,handles.ray_number);
handles.weighting_factor=weighting_factor;



p=awgn(reshape(handles.projection_data,1,[])',handles.noise_level);
f_temp=zeros(handles.image_size^2,1);
f_final=zeros(handles.image_size^2,1);


steps = handles.iteration_number;
for ii=1:handles.iteration_number
    
 [result_image,f_final,f_temp]=Func_ART(p,weighting_factor,handles.relaxation_factor,f_final,f_temp,handles.image_size);
 
 
 if getappdata(h,'canceling')
        break
    end

    waitbar(ii/handles.iteration_number,h,[num2str(ii) ' iteration done'])
 
end

delete(h) 


handles.result_image=result_image;
guidata(hObject,handles);



axes(handles.result)
imagesc(result_image)

title(['M=', num2str(handles.image_size),', N=',num2str(handles.iteration_number),', \lambda=',num2str(handles.relaxation_factor),',  SNR=',num2str(handles.noise_level), ' dB ' ,', w-' handles.w_name ]) 






% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)

result_image=handles.result_image;
uisave('result_image','result_image')
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





% Hints: get(hObject,'String') returns contents of snr as text
%        str2double(get(hObject,'String')) returns contents of snr as a double


% --- Executes during object creation, after setting all properties.
function snr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to snr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in popup.
function popup_Callback(hObject, eventdata, handles)
result_image=handles.result_image;

figure
imagesc(result_image);
title(['M=', num2str(handles.image_size),', N=',num2str(handles.iteration_number),', \lambda=',num2str(handles.relaxation_factor),',  SNR=',num2str(handles.noise_level), ' dB ' ,', w-' handles.w_name ])

% hObject    handle to popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
set(handles.iteration,'string','');
set(handles.lambda,'string','');
set(handles.snr,'string','');
set(handles.alpha,'value',1);
set(handles.size,'value',1);

cla(handles.result,'reset');
cla(handles.sinogram,'reset');

set(handles.size,'Enable','off');
set(handles.iteration,'Enable','off');
set(handles.lambda,'Enable','off');
set(handles.alpha,'Enable','off');
set(handles.snr,'Enable','off');
set(handles.reconstruct,'Enable','off');

% clear all;

% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
% function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3
