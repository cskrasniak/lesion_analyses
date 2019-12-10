
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
                     
s0 = daq.createSession('ni');
AI2 = addAnalogInputChannel(s0,'dev1','ai2', 'Voltage');
s0.Rate = 100000;
s0.IsContinuous = true;
%s0.DurationInSeconds = 60;
s0.IsNotifyWhenDataAvailableExceedsAuto = true;


%AI0=addAnalogInputChannel(s0,'ai0', 'Voltage');
%AI1=addAnalogInputChannel(s0,'dev1','ai1', 'Voltage');
% s1.Rate = 8000;
% s1.IsContinuous = true;
% s1.startBackground
global s2
s2 = daq.createSession('ni');
AO0=addAnalogOutputChannel(s2,'dev1','ao0', 'Voltage');
AO1=addAnalogOutputChannel(s2,'dev1','ao1', 'Voltage');
A02 = addAnalogOutputChannel(s2,'dev2','ao1', 'Voltage');
s2.Rate = 8000;

risetime = 0.01; % rise time in seconds
durxy = .05;% Duration time for signal in seconds
xConv = 0.175;% Conversion rate from mm to volt x axis
yConv = 0.1675;% Conversion rate from mm to volt y axis
XY_list = [0,0];

%fid1 = fopen('data.bin','w');
formatOut = 'yyyy-mm-dd';
date = datestr(now,formatOut);
dataPath = string(mousedir)+'\'+date;
mkdir(dataPath);cd(dataPath)
saveName = mouseName+"_"+date+"_1";
%% laser stimulation specs
dt = 1/s2.Rate;%seconds
stopTime2 = 2; %downward amplitude ramp period/length of trial
stopTime1 = stopTime2-.1; %seconds
t1 = 0:dt:stopTime1-dt;
t2 = stopTime1:dt:stopTime2-dt;
ampmax = 5;% needs to be set to 1
ampRamp = linspace(0,ampmax,length(t2));
amp=repmat(ampmax,1,length(ampRamp))-ampRamp;

lo1 = ampmax*sin(2*pi*t1*40)+ampmax; % front number is amplitude, 40 is 40hz stim, last is to make it all positive
lo2 = amp.*sin(2*pi*t2*40)+amp;
laserOutput = [lo1,lo2];
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

%fclose(fid1);
XY_list = XY_list(2:end,:);
input = inputdlg("was the laser on? yes=1, no = 0","laser on?");
XY_list(:,3) = repmat(str2double(input{1}),size(XY_list,1),1);
XY_list(:,2) = XY_list(:,2).*-1; %make the y values negative b/c pos/neg is switched for this mirror
if exist(saveName,'file') %save the file
    save(saveName,'XY_list')
else
    num(1) = 0;
    filelist=dir('*.mat');%if one experiment has already been done on this mouse today, save it under the next number
    for i= 1:length(filelist)
        num(i) = str2num(filelist(i).name(end-4));
    end

    save(mouseName+"_"+date+"_"+string(max(num)+1),'XY_list');
end
clear num





function [X,Y] = getLocation(s)
[xin, yin] = inputSingleScan(s);
end
function logData(src,event,fid)
     fwrite(fid,event.Data, 'double');
%     save('data', event.Data)
end

function [X,Y] = getDestination(XY_list)
    % destListx = [repmat(-3,5,1);repmat(-2,7,1);repmat(-1,7,1);repmat(0,9,1);repmat(-1,9,1);repmat(-2,9,1);repmat(-3,9,1)];%these are all the same locations used in guo et al, 2014
    % destListy = [[-2:2]';[-3:3]';[-3:3]';[-4:4]';[-4:4]';[-4:4]';[-4:4]' ];
    destListx = [[-2.5:1:2.5],[-2.5:1:2.5],[-3.5:1:3.5],[-3.5:1:3.5],[-3.5:1:3.5],[-3.5:1:3.5],[-3.5:1:3.5],[-3.5:1:3.5],[-3.5:1:3.5],[-1.5:1:1.5],2.5,-2.5]';%these are based off my surgeries
    destListy = [repmat(-3,6,1);repmat(-2,6,1);repmat(-1,8,1);repmat(0,8,1);repmat(1,8,1);repmat(2,8,1);repmat(3,8,1);repmat(4,8,1);repmat(5,8,1);repmat(6,4,1);repmat(8,2,1)];% the last two numbers are on the headplate as negative controls, %the y mirror is opposite so need to flip it
    destList = [destListx,destListy];
    idx=randsample([1:size(destList,1)],1,true);
    X = destList(idx,1); Y = destList(idx,2);
    if X == XY_list(end,1) && Y == XY_list(end,2)
        [X,Y] = getDestination(XY_list);
    end
end