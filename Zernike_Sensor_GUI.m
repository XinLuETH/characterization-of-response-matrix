function varargout = Zernike_Sensor_GUI(varargin)
% ZERNIKE_SENSOR_GUI MATLAB code for Zernike_Sensor_GUI.fig
%      ZERNIKE_SENSOR_GUI, by itself, creates a new ZERNIKE_SENSOR_GUI or raises the existing
%      singleton*.
%
%      H = ZERNIKE_SENSOR_GUI returns the handle to a new ZERNIKE_SENSOR_GUI or the handle to
%      the existing singleton*.
%
%      ZERNIKE_SENSOR_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ZERNIKE_SENSOR_GUI.M with the given input arguments.
%
%      ZERNIKE_SENSOR_GUI('Property','Value',...) creates a new ZERNIKE_SENSOR_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Zernike_Sensor_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Zernike_Sensor_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Zernike_Sensor_GUI

% Last Modified by GUIDE v2.5 22-Feb-2017 11:00:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Zernike_Sensor_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Zernike_Sensor_GUI_OutputFcn, ...
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


% --- Executes just before Zernike_Sensor_GUI is made visible.
function Zernike_Sensor_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Zernike_Sensor_GUI (see VARARGIN)

% Choose default command line output for Zernike_Sensor_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Zernike_Sensor_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Zernike_Sensor_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in read_1fits.
function read_1fits_Callback(hObject, eventdata, handles)
% hObject    handle to read_1fits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.fits','Select the fits file');
if FileName ~= 0
   data0 = fitsread(strcat(PathName, FileName)); 
   axes(handles.Image);
   imagesc(data0(:,:,1))%
   colormap Jet
   set(gca,'YDir','normal')
   hold on
end



function X_Callback(hObject, eventdata, handles)
% hObject    handle to X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of X as text
%        str2double(get(hObject,'String')) returns contents of X as a double


% --- Executes during object creation, after setting all properties.
function X_CreateFcn(hObject, eventdata, handles)
% hObject    handle to X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Y_Callback(hObject, eventdata, handles)
% hObject    handle to Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Y as text
%        str2double(get(hObject,'String')) returns contents of Y as a double


% --- Executes during object creation, after setting all properties.
function Y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R_Callback(hObject, eventdata, handles)
% hObject    handle to R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R as text
%        str2double(get(hObject,'String')) returns contents of R as a double


% --- Executes during object creation, after setting all properties.
function R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_circle.
function plot_circle_Callback(hObject, eventdata, handles)
% hObject    handle to plot_circle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X = str2double(get(handles.X, 'String'));
Y = str2double(get(handles.Y, 'String'));
radius = str2double(get(handles.R, 'String'));
center = [X,Y];
axes(handles.Image);
plotH = viscircles(center,radius);
pause(5) 
delete(plotH)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
% read fits files
[FileName,PathName] = uigetfile('*.fits','Select the fits file of 0*Pi/2');
if FileName ~= 0
   data0 = fitsread(strcat(PathName, FileName)); 
end

[FileName,PathName] = uigetfile('*.fits','Select the fits file of 1*Pi/2');
if FileName ~= 0
   data1 = fitsread(strcat(PathName, FileName)); 
end

[FileName,PathName] = uigetfile('*.fits','Select the fits file of 2*Pi/2');
if FileName ~= 0
   data2 = fitsread(strcat(PathName, FileName)); 
end

[FileName,PathName] = uigetfile('*.fits','Select the fits file of 3*Pi/2');
if FileName ~= 0
   data3 = fitsread(strcat(PathName, FileName)); 
end

X = str2double(get(handles.X, 'String'));
Y = str2double(get(handles.Y, 'String'));
radius = str2double(get(handles.R, 'String'));

% reduce the data size
DATA0 = data0(Y-radius-5:Y+radius+5, X-radius-5:X+radius+5,:);
DATA1 = data1(Y-radius-5:Y+radius+5, X-radius-5:X+radius+5,:);
DATA2 = data2(Y-radius-5:Y+radius+5, X-radius-5:X+radius+5,:);
DATA3 = data3(Y-radius-5:Y+radius+5, X-radius-5:X+radius+5,:);
num_points = size(DATA0,1)
% get the average of DATA0 in the ROI and normalize the axes
[x, y]=meshgrid(linspace((-radius-5)/radius,(radius+5)/radius, num_points),linspace((-radius-5)/radius,(radius+5)/radius, num_points));
[qi,ri] = cart2pol(x,y);
IOI = ri<=1;
% the cernter of the circle changes after resize the data
% X = radius + 5 + 1;
% Y = radius + 5 + 1;
% This assume the circle falls *entirely* inside the image
Circle = (x-0).^2+(y-0).^2 <= 1^2;
DATA0(~Circle) = 0;
DATA1(~Circle) = 0;
DATA2(~Circle) = 0;
DATA3(~Circle) = 0;

%% remove the patters 
% set data outside the circle to zero
% for i = 1:size(DATA0,3)
%     f0 = fftshift(fft2(DATA0(:,:,i)));
%     f1 = fftshift(fft2(DATA1(:,:,i)));
%     f2 = fftshift(fft2(DATA2(:,:,i)));
%     f3 = fftshift(fft2(DATA3(:,:,i)));
% % used to plot the image
%     fLog0 = log(1 + abs(f0));
%     fLog1 = log(1 + abs(f1));
%     fLog2 = log(1 + abs(f2));
%     fLog3 = log(1 + abs(f3));
% % filter by a range based on fLog
%     filter0 = fLog0 < .7*max(fLog0(:));
%     filter1 = fLog1 < .7*max(fLog1(:));
%     filter2 = fLog2 < .7*max(fLog2(:));
%     filter3 = fLog3 < .7*max(fLog3(:));
%     filter0(54:78,54:78) = 1;
%     filter1(54:78,54:78) = 1;
%     filter2(54:78,54:78) = 1;
%     filter3(54:78,54:78) = 1;
%     DATA0(:,:,i) = abs(ifft2(f0.*filter0));
%     DATA1(:,:,i) = abs(ifft2(f1.*filter1));
%     DATA2(:,:,i) = abs(ifft2(f2.*filter2));
%     DATA3(:,:,i) = abs(ifft2(f3.*filter3));
% end

%% Zernike sensing
contents = cellstr(get(hObject,'String'));
if strcmp(contents{get(hObject,'Value')}, 'Phase shifting')
    I0 = (DATA0+DATA1+DATA2+DATA3)/4;
    phase = (DATA1-DATA3)./I0;
    phase = phase*633/2/pi;
    phase(~Circle) = 0;
    
%     a = zernike_coeffs(phase(:,:,1), 15);
    
elseif strcmp(contents{get(hObject,'Value')}, 'Standard method')
    b = 0.5;
    % normalize other datas to data0
    for i=1:size(DATA0,3)
        A = DATA0(:,:,i); % choose one image
        val = mean(A(Circle)); % mean
        DATA1(:,:,i) = DATA1(:,:,i)/val;
        DATA2(:,:,i) = DATA2(:,:,i)/val;
        DATA3(:,:,i) = DATA3(:,:,i)/val;
    end
    phase = (-1+sqrt(3-2*b-(1-DATA1)/b))/pi;
    phase = phase*633/2/pi;
    phase(~Circle) = 0;
%     a = zernike_coeffs(phase(:,:,1), 15);
end
wavefront = phase(:,:,1);
axes(handles.Wavefront);
imagesc(phase(:,:,1))
set(gca,'YDir','normal')
colormap Jet

% get the first 15 zernike terms
n = [0  1 1  2 2 2  3  3 3 3  4  4 4 4 4];
m = [0 -1 1 -2 0 2 -3 -1 1 3 -4 -2 0 2 4];

Z = zernfun(n,m,ri(IOI),qi(IOI),'norm');
% decompose the wavefront
a = Z\wavefront(IOI)

% send the coefficient to each "term"

textLabel = sprintf('%s', num2str(a(1)));
set(handles.term1, 'String', textLabel);
textLabel = sprintf('%s', num2str(a(2)));
set(handles.term2, 'String', textLabel);
textLabel = sprintf('%s', num2str(a(3)));
set(handles.term3, 'String', textLabel);
textLabel = sprintf('%s', num2str(a(4)));
set(handles.term4, 'String', textLabel);
textLabel = sprintf('%s', num2str(a(5)));
set(handles.term5, 'String', textLabel);
textLabel = sprintf('%s', num2str(a(6)));
set(handles.term6, 'String', textLabel);
textLabel = sprintf('%s', num2str(a(7)));
set(handles.term7, 'String', textLabel);
textLabel = sprintf('%s', num2str(a(8)));
set(handles.term8, 'String', textLabel);
textLabel = sprintf('%s', num2str(a(9)));
set(handles.term9, 'String', textLabel);
textLabel = sprintf('%s', num2str(a(10)));
set(handles.term10, 'String', textLabel);
textLabel = sprintf('%s', num2str(a(11)));
set(handles.term11, 'String', textLabel);
textLabel = sprintf('%s', num2str(a(12)));
set(handles.term12, 'String', textLabel);
textLabel = sprintf('%s', num2str(a(13)));
set(handles.term13, 'String', textLabel);
textLabel = sprintf('%s', num2str(a(14)));
set(handles.term14, 'String', textLabel);
textLabel = sprintf('%s', num2str(a(15)));
set(handles.term15, 'String', textLabel);




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



function term1_Callback(hObject, eventdata, handles)
% hObject    handle to term1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of term1 as text
%        str2double(get(hObject,'String')) returns contents of term1 as a double


% --- Executes during object creation, after setting all properties.
function term1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to term1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function term2_Callback(hObject, eventdata, handles)
% hObject    handle to term2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of term2 as text
%        str2double(get(hObject,'String')) returns contents of term2 as a double


% --- Executes during object creation, after setting all properties.
function term2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to term2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function term3_Callback(hObject, eventdata, handles)
% hObject    handle to term3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of term3 as text
%        str2double(get(hObject,'String')) returns contents of term3 as a double


% --- Executes during object creation, after setting all properties.
function term3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to term3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function term4_Callback(hObject, eventdata, handles)
% hObject    handle to term4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of term4 as text
%        str2double(get(hObject,'String')) returns contents of term4 as a double


% --- Executes during object creation, after setting all properties.
function term4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to term4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function term5_Callback(hObject, eventdata, handles)
% hObject    handle to term5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of term5 as text
%        str2double(get(hObject,'String')) returns contents of term5 as a double


% --- Executes during object creation, after setting all properties.
function term5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to term5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function term6_Callback(hObject, eventdata, handles)
% hObject    handle to term6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of term6 as text
%        str2double(get(hObject,'String')) returns contents of term6 as a double


% --- Executes during object creation, after setting all properties.
function term6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to term6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function term7_Callback(hObject, eventdata, handles)
% hObject    handle to term7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of term7 as text
%        str2double(get(hObject,'String')) returns contents of term7 as a double


% --- Executes during object creation, after setting all properties.
function term7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to term7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function term8_Callback(hObject, eventdata, handles)
% hObject    handle to term8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of term8 as text
%        str2double(get(hObject,'String')) returns contents of term8 as a double


% --- Executes during object creation, after setting all properties.
function term8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to term8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function term9_Callback(hObject, eventdata, handles)
% hObject    handle to term9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of term9 as text
%        str2double(get(hObject,'String')) returns contents of term9 as a double


% --- Executes during object creation, after setting all properties.
function term9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to term9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function term10_Callback(hObject, eventdata, handles)
% hObject    handle to term10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of term10 as text
%        str2double(get(hObject,'String')) returns contents of term10 as a double


% --- Executes during object creation, after setting all properties.
function term10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to term10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function term11_Callback(hObject, eventdata, handles)
% hObject    handle to term11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of term11 as text
%        str2double(get(hObject,'String')) returns contents of term11 as a double


% --- Executes during object creation, after setting all properties.
function term11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to term11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function term12_Callback(hObject, eventdata, handles)
% hObject    handle to term12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of term12 as text
%        str2double(get(hObject,'String')) returns contents of term12 as a double


% --- Executes during object creation, after setting all properties.
function term12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to term12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function term13_Callback(hObject, eventdata, handles)
% hObject    handle to term13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of term13 as text
%        str2double(get(hObject,'String')) returns contents of term13 as a double


% --- Executes during object creation, after setting all properties.
function term13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to term13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function term14_Callback(hObject, eventdata, handles)
% hObject    handle to term14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of term14 as text
%        str2double(get(hObject,'String')) returns contents of term14 as a double


% --- Executes during object creation, after setting all properties.
function term14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to term14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function term15_Callback(hObject, eventdata, handles)
% hObject    handle to term15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of term15 as text
%        str2double(get(hObject,'String')) returns contents of term15 as a double


% --- Executes during object creation, after setting all properties.
function term15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to term15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
