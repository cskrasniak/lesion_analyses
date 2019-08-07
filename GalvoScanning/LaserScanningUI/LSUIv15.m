function varargout = LSUIv15(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Copyrights 2018, Hemanth Mohan, Cold Spring Harbor Labs. email: mohan (at) cshl (dot) edu %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LSUIv15_OpeningFcn, ...
                   'gui_OutputFcn',  @LSUIv15_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before LSUIv15 is made visible.
function LSUIv15_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
set(handles.Vstatus,'String','Ready')
%% setup Galvos
handles.s = daq.createSession('ni');
addAnalogOutputChannel(handles.s,'Dev1',1,'Voltage');
addAnalogOutputChannel(handles.s,'Dev1',0,'Voltage');
handles.s.Rate = 1000;
handles.xyval = [str2double(get(handles.x_volt,'String')) str2double(get(handles.y_volt,'String'))];
%% setup Laser Pulse
% handles.s2 = daq.createSession('ni');
% handles.s2.Rate = 1000;
% addCounterOutputChannel(handles.s2,'dev1', 0, 'PulseGeneration');
%% setup reading trigger
handles.s3 = daq.createSession('ni');
addAnalogInputChannel(handles.s3,'Dev1',1,'Voltage');
addAnalogInputChannel(handles.s3,'Dev1',0,'Voltage');
handles.DataTable = cell2table(cell(0,8),'VariableNames',{'VidID', ...
    'Trigger','Xloc','Yloc','Freq','IDelay','Duration','DCycle'});

%% Video Preview
% imaqreset;
% handles.vid = videoinput('gige', 1, 'Mono8');
% src = getselectedsource(handles.vid);
% src.AcquisitionFrameRateEnable = 'True';
% src.AcquisitionFrameRate = 11;
% vidRes = handles.vid.VideoResolution;
% nBands = handles.vid.NumberOfBands;
% handles.hImage = image(handles.VideoDisp , zeros(vidRes(2), vidRes(1), nBands));
% handles.output = hObject;
% Update handles structure %ask Hemanth for help here
guidata(hObject, handles);

% UIWAIT makes LSUIv15 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LSUIv15_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function getvolt_Callback(hObject, eventdata, handles)
handles.oldVoltx = str2double(get(handles.x_volt,'String'));
handles.oldVolty = str2double(get(handles.y_volt,'String'));
[xval yval] = ginput(1);
set(handles.x_volt,'String',xval);
set(handles.y_volt,'String',yval);

guidata(hObject,handles)


% --- Executes on button press in getvolt.

% hObject    handle to getvolt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in StartV.
function StartV_Callback(hObject, eventdata, handles)

% hObject    handle to StartV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.xyval = [handles.xyval; str2double( get(handles.x_volt,'String')) str2double( get(handles.y_volt,'String'))];
set(handles.Vstatus,'String','Wait')

%% Modeling input signal %%
risetime = 0.1; % rise time in seconds
durxy = 1;% Duration time for signal in seconds
xConv = .35/2;% Conversion rate from mm to volt x axis  .4
Ampx = str2double( get(handles.x_volt,'String')) * xConv;
durx_temp = durxy; %%% Duration in Seconds
durx = 0:1/handles.s.Rate:durx_temp;
outx_1 = Ampx * ones(1,size(durx,2)); 
rsx1 = find(durx<risetime);
risetimesize = size(rsx1,2);
rsx2 = linspace(handles.xyval(end-1,1),Ampx,size(rsx1,2)); %%%% set values between [ ] to manage sigmoid risetime
outx_1(1:size(rsx2,2)) = rsx2; 
outputx = outx_1';

yConv = .335/2;% Conversion rate from mm to volt y axis  .39
Ampy = str2double( get(handles.y_volt,'String')) * yConv;
dury_temp = durxy; %%% Duration in Seconds
dury = 0:1/handles.s.Rate:dury_temp;
outy_1 = Ampy * ones(1,size(dury,2)); 
rsy1 = find(dury<risetime);
rsy2 = linspace(handles.xyval(end-1,2),Ampy,size(rsy1,2)); %%%% set values between [ ] to manage sigmoid risetime
outy_1(1:size(rsy2,2)) = rsy2; 
outputy = outy_1';


%%tan Output data to daq %%
queueOutputData(handles.s,[outputx outputy]);
tic
handles.s.startForeground;
OutputDuration = toc;
outputx = [];
outputy = [];

disp(['The output duration is: ' mat2str(OutputDuration)]);
set(handles.Vstatus,'String','Ready')
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1



function y_volt_Callback(hObject, eventdata, handles)
% hObject    handle to y_volt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_volt as text
%        str2double(get(hObject,'String')) returns contents of y_volt as a double
yvolt = str2double(get(hObject,'String'));
if yvolt>10 || yvolt<-10
    set(handles.y_volt,'String','Error')
end

% --- Executes during object creation, after setting all properties.
function y_volt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_volt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function x_volt_Callback(hObject, eventdata, handles)
% hObject    handle to x_volt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_volt as text
%        str2double(get(hObject,'String')) returns contents of x_volt as a double
xvolt = str2double(get(hObject,'String'));
if xvolt>10 || xvolt<-10
    set(handles.x_volt,'String','Error')
end

% --- Executes during object creation, after setting all properties.

function x_volt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_volt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes during object creation, after setting all properties.
function lasernav_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lasernav (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate lasernav



function Vstatus_Callback(hObject, eventdata, handles)
% hObject    handle to Vstatus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vstatus as text
%        str2double(get(hObject,'String')) returns contents of Vstatus as a double


% --- Executes during object creation, after setting all properties.
function Vstatus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vstatus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in StartVid.
function StartVid_Callback(hObject, eventdata, handles)
% hObject    handle to StartVid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

preview(handles.vid, handles.hImage);
handles.VideoDisp.PlotBoxAspectRatioMode = 'manual';
handles.VideoDisp.PlotBoxAspectRatio = [1312/1082, 1,1];  %%%% 1312/1082
set(handles.VidStatus,'String','Video Running')

% --- Executes during object creation, after setting all properties.
function VideoDisp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VideoDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate VideoDisp


% --- Executes on button press in stopvid.
function stopvid_Callback(hObject, eventdata, handles)
% hObject    handle to stopvid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stoppreview(handles.vid);
set(handles.VidStatus,'String','Video Stopped')


% --- Executes on button press in SnapShot.
function SnapShot_Callback(hObject, eventdata, handles)
% hObject    handle to SnapShot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
frame = getsnapshot(handles.vid);
image(handles.lasernav,frame);
handles.lasernav.PlotBoxAspectRatioMode = 'manual';
handles.lasernav.PlotBoxAspectRatio = [1312/1082, 1,1];  %%% 8 um per pixel
handles.lasernav.XTick =[];
handles.lasernav.YTick =[];



function Freq_Callback(hObject, eventdata, handles)
% hObject    handle to Freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Freq as text
%        str2double(get(hObject,'String')) returns contents of Freq as a double


% --- Executes during object creation, after setting all properties.
function Freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Duty_Callback(hObject, eventdata, handles)
% hObject    handle to Duty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Duty as text
%        str2double(get(hObject,'String')) returns contents of Duty as a double


% --- Executes during object creation, after setting all properties.
function Duty_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Duty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function IDelay_Callback(hObject, eventdata, handles)
% hObject    handle to IDelay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IDelay as text
%        str2double(get(hObject,'String')) returns contents of IDelay as a double


% --- Executes during object creation, after setting all properties.
function IDelay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IDelay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LUpdate.
function LUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to LUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.s2.Channels.Frequency = str2double( get(handles.Freq,'String'));
handles.s2.Channels.DutyCycle = str2double( get(handles.Duty,'String'));
handles.s2.Channels.InitialDelay = str2double( get(handles.IDelay,'String'));
handles.s2.DurationInSeconds = str2double( get(handles.IDelay,'String')) + str2double( get(handles.Duration,'String'));

Pt = 0 : 1/handles.s2.Rate : handles.s2.DurationInSeconds;        
PulsPer = 0 : 1/handles.s2.Channels.Frequency : handles.s2.DurationInSeconds; 
PulWidth = handles.s2.Channels.DutyCycle/handles.s2.Channels.Frequency;
Signal = pulstran(Pt,PulsPer,'rectpuls', PulWidth);
area(handles.axes1,Pt, Signal);
axis(handles.axes1,[0 1 0 1])

% --- Executes on button press in LOn.
function LOn_Callback(hObject, eventdata, handles)
% hObject    handle to LOn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.s2.IsContinuous = 1;
startBackground(handles.s2);


% --- Executes on button press in LOff.
function LOff_Callback(hObject, eventdata, handles)
% hObject    handle to LOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stop(handles.s2)
handles.s2.IsContinuous = 0;



function Duration_Callback(hObject, eventdata, handles)
% hObject    handle to Duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Duration as text
%        str2double(get(hObject,'String')) returns contents of Duration as a double


% --- Executes during object creation, after setting all properties.
function Duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function VidStatus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VidStatus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in trigger.
function trigger_Callback(hObject, eventdata, handles)
% hObject    handle to trigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

startBackground(handles.s2);


% --- Executes on button press in readtrigger.
function readtrigger_Callback(hObject, eventdata, handles)
% hObject    handle to readtrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of readtrigger
while get(hObject,'Value')
    Sig = handles.s3.inputSingleScan;
    if Sig(1) > 4
        startBackground(handles.s2);
        CurID = str2num(get(handles.CurrentID,'string'));
        handles.DataTable.VidID(1,1) = CurID;
        handles.DataTable.Trigger(1,1) = 1;
        handles.DataTable.Xloc(1,1) = str2num(get(handles.x_volt,'string'));
        handles.DataTable.Yloc(1,1) = str2num(get(handles.y_volt,'string'));
        handles.DataTable.Freq(1,1) = str2num(get(handles.Freq,'string'));
        handles.DataTable.IDelay(1,1) = str2num(get(handles.IDelay,'string'));
        handles.DataTable.Duration(1,1) = str2num(get(handles.Duration,'string'));
        handles.DataTable.DCycle(1,1) = str2num(get(handles.Duty,'string'));
        StimulusTable = handles.DataTable;
        if exist(fullfile([get(handles.FilePath,'string')], ['StimTable.mat'])) == 2
            STable = load(fullfile([get(handles.FilePath,'string')], ['StimTable.mat']));
            StimulusTable = [STable.StimulusTable; StimulusTable];
            save(fullfile([get(handles.FilePath,'string')], ['StimTable.mat']),'StimulusTable');
        else
            save(fullfile([get(handles.FilePath,'string')], ['StimTable.mat']),'StimulusTable');
        end
        pause(handles.s2.DurationInSeconds - (handles.s2.DurationInSeconds/2))
        CurID = CurID+1;
        set(handles.CurrentID,'string',CurID);
    elseif Sig(2) > 4
        CurID = str2num(get(handles.CurrentID,'string'));
        handles.DataTable.VidID(1,1) = CurID;
        handles.DataTable.Trigger(1,1) = 0;
        handles.DataTable.Xloc(1,1) = str2num(get(handles.x_volt,'string'));
        handles.DataTable.Yloc(1,1) = str2num(get(handles.y_volt,'string'));
        handles.DataTable.Freq(1,1) = str2num(get(handles.Freq,'string'));
        handles.DataTable.IDelay(1,1) = str2num(get(handles.IDelay,'string'));
        handles.DataTable.Duration(1,1) = str2num(get(handles.Duration,'string'));
        handles.DataTable.DCycle(1,1) = str2num(get(handles.Duty,'string'));
        StimulusTable = handles.DataTable;
        if exist(fullfile([get(handles.FilePath,'string')], ['StimTable.mat'])) == 2
            STable = load(fullfile([get(handles.FilePath,'string')], ['StimTable.mat']));
            StimulusTable = [STable.StimulusTable; StimulusTable];
            save(fullfile([get(handles.FilePath,'string')], ['StimTable.mat']),'StimulusTable');
        else
            save(fullfile([get(handles.FilePath,'string')], ['StimTable.mat']),'StimulusTable');
        end
        pause(handles.s2.DurationInSeconds - (handles.s2.DurationInSeconds/2))
        CurID = CurID+1;
        set(handles.CurrentID,'string',CurID);
    end
    drawnow();
end



function FilePath_Callback(hObject, eventdata, handles)
% hObject    handle to FilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FilePath as text
%        str2double(get(hObject,'String')) returns contents of FilePath as a double


% --- Executes during object creation, after setting all properties.
function FilePath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EnterID_Callback(hObject, eventdata, handles)
% hObject    handle to EnterID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EnterID as text
%        str2double(get(hObject,'String')) returns contents of EnterID as a double


% --- Executes during object creation, after setting all properties.
function EnterID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EnterID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in updateID.
function updateID_Callback(hObject, eventdata, handles)
% hObject    handle to updateID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.CurrentID,'string',get(handles.EnterID,'string'));

% --- Executes during object creation, after setting all properties.
function CurrentID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrentID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
