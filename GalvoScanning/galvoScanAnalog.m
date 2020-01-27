%%% Script for running galvo scanning in the IBL task, ran in unison with
%%% biasedScanningChoiceWorld in pyBpod. first execute this script, then
%%% start the task in pybpod. When ending the session, stop this script by
%%% pressing the big "end now" button then stop the pybpod task. Written
%%% for use with two NI-USB-6211 DAQ boards. Written by Christopher S
%%% Krasniak, Cold Spring Harbor Laboratory/International Brain Lab, August
%%% 2019.



%% End trial button
mousedir = uigetdir('C:\Users\IBLuser\Documents\laserPostitionData',"What's the mouse's name?");
mouseDirSplit = strsplit(mousedir,'\');
mouseName = mouseDirSplit{end};
ButtonHandle = uicontrol('Style', 'PushButton', ...
                         'String', 'End Now', ...
                         'Position', [10,10,500,400],...
                         'FontSize', 48,...
                         'BackgroundColor', 'red',...
                         'ForegroundColor','white',...
                         'Callback', 'delete(gcbf)');
%% Setting up DAQ boards                     
s0 = daq.createSession('ni');
AI2 = addAnalogInputChannel(s0,'dev1','ai2', 'Voltage');
s0.Rate = 100000;
s0.IsContinuous = true;
s0.IsNotifyWhenDataAvailableExceedsAuto = true;
s0.NotifyWhenDataAvailableExceeds = 1000;

AI0=addAnalogInputChannel(s0,'dev1','ai0', 'Voltage');
AI1=addAnalogInputChannel(s0,'dev1','ai1', 'Voltage');

global s2
s2 = daq.createSession('ni');
AO0=addAnalogOutputChannel(s2,'dev1','ao0', 'Voltage');
AO1=addAnalogOutputChannel(s2,'dev1','ao1', 'Voltage');
A02 = addAnalogOutputChannel(s2,'dev2','ao1', 'Voltage');
s2.Rate = 8000;

%% saving Parameters
formatOut = 'yyyy-mm-dd';
date = datestr(now,formatOut);
dataPath = string(mousedir)+'\'+date;
mkdir(dataPath);cd(dataPath)
saveName = mouseName+"_"+date+"_1";

%% laser stimulation specs
xConv = 0.175;% Conversion rate from mm to volt x axis
yConv = 0.1675;% Conversion rate from mm to volt y axis
XY_list = [0,0];
dt = 1/s2.Rate;%seconds
stopTime2 = 2; %downward amplitude ramp period/length of trial
stopTime1 = stopTime2-.1; %seconds
t1 = 0:dt:stopTime1-dt;
t2 = stopTime1:dt:stopTime2-dt; % time for the ramp down
ampmax = 5;% needs to be set to 1
ampRamp = linspace(0,ampmax,length(t2)); 
amp=repmat(ampmax,1,length(ampRamp))-ampRamp; % the final amplitudes for the ramp down

lo1 = ampmax*sin(2*pi*t1*40)+ampmax; % front number is amplitude, 40 is 40hz stim, last is to make it all positive
lo2 = amp.*sin(2*pi*t2*40)+amp; %creating the downward ramp of the sine amp
laserOutput = [lo1,lo2]; % put together normal output and ramp
%% Main loop
while true
    tic

  if ~ishandle(ButtonHandle)
      release(s0)
    disp('experiment stopped by user');
    break;
  end

[X_dest,Y_dest] = getDestination(XY_list);
Ampx = X_dest*xConv;
outputx = Ampx';


Ampy = Y_dest * yConv;
outputy = Ampy';
 %%%%%%%%%%%%%%%%%%% NEED Y THEN X FOR OUTPUT
queueOutputData(s2,[repmat(outputy,length(laserOutput),1) repmat(outputx,length(laserOutput),1) laserOutput'])

newTrialListener = addlistener(s0,'DataAvailable', @newTrialCheck); %Add listener to check if there is a new trial aka if the laser should move 

s0.startBackground
toc
try
    fprintf('Moving to %.2f,%.2f \n ',X_dest,Y_dest*-1)
    s0.wait(61)%hold all operations for 61 seconds, (1s more than max trial length), if new trial not triggered by then, this trial is repeated
    XY_list = [XY_list; [X_dest,Y_dest]]; % the Y is backwards so need to negate it
catch
    s2.startForeground 
    s0.stop()
    delete(newTrialListener); %delete(mirrorPosListener);
    disp("something went wrong, no trigger detected")
end
delete(newTrialListener); %delete(mirrorPosListener);

end
%% reset laser position to 0
disp("Moving to 0,0")
queueOutputData(s2,zeros(100,3))
s2.startForeground()
s2.release()
s0.release()

%% Saving data 
XY_list = XY_list(2:end,:);
input = inputdlg("was the laser on? yes=1, no = 0","laser on?");
XY_list(:,3) = repmat(str2double(input{1}),size(XY_list,1),1);
XY_list(:,2) = XY_list(:,2).*-1; %make the y values negative b/c pos/neg is switched for this mirror
if exist(saveName,'file') %save the file
    save(saveName,'XY_list');
    disp("Data Saved")
else
    num(1) = 0;
    filelist=dir('*.mat');%if one experiment has already been done on this mouse today, save it under the next number
    for i= 1:length(filelist)
        num(i) = str2num(filelist(i).name(end-4));
    end

    save(mouseName+"_"+date+"_"+string(max(num)+1),'XY_list');
    %writeNPY('XY_list',mouseName+"_"+date+"_"+string(max(num)+1)+".npy");
    disp("Data Saved")
end
clear num

%% Helper functions

function [X,Y] = getDestination(XY_list)

    destListx = [[-1.5:1:1.5],[-2.5:1:2.5],[-3.5:1:3.5],[-3.5:1:3.5],[-4.5:1:4.5],[-4.5:1:4.5],[-4.5:1:4.5],[-3.5:1:3.5],2.5,-2.5,2.5,-2.5]';%these are based off my surgeries
    % the last two numbers are on the headplate as negative controls, 
    %the y mirror is backwards so need to flip it
    % subtract .5 to move it back .5cm, off OB
    destListy = [repmat(-3,4,1);repmat(-2,6,1);repmat(-1,8,1);repmat(0,8,1);repmat(1,10,1);repmat(2,10,1);repmat(3,10,1);repmat(4,8,1);repmat(8,4,1)] -.5;
    destList = [destListx,destListy];
    idx=randsample([1:size(destList,1)],1,true);
    X = destList(idx,1); Y = destList(idx,2);
    if X == XY_list(end,1) && Y == XY_list(end,2)
        [X,Y] = getDestination(XY_list);
    end
end