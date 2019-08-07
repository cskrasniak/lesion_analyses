
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
%s0.IsNotifyWhenDataAvailableExceedsAuto = true;
s0.NotifyWhenDataAvailableExceeds = 1000;

AI0=addAnalogInputChannel(s0,'dev1','ai0', 'Voltage');
AI1=addAnalogInputChannel(s0,'dev1','ai1', 'Voltage');
global s2
s2 = daq.createSession('ni');
AO0=addAnalogOutputChannel(s2,'dev1','ao0', 'Voltage');
AO1=addAnalogOutputChannel(s2,'dev1','ao1', 'Voltage');
s2.Rate = 80000;

risetime = 0.01; % rise time in seconds
durxy = .05;% Duration time for signal in seconds
xConv = 0.175;% Conversion rate from mm to volt x axis
yConv = 0.1675;% Conversion rate from mm to volt y axis

%% Save parameters
formatOut = 'yyyy-mm-dd';
date = datestr(now,formatOut);
dataPath = string(mousedir)+'\'+date;
mkdir(dataPath);cd(dataPath)
saveName = mouseName+"_"+date+"_1";

%% generate a list to use for stimulation 
visCtx_x = [1.5,1.5,2.5,2.5,-1.5,-1.5,-2.5,-2.5];
visCtx_y = [6,5,5,6,6,5,5,6];
offTarget_x = [1.5,1.5,2.5,2.5,-1.5,-1.5,-2.5,-2.5];%list of false stimulation sites, must be the same size as target_x/y
offTarget_y = [10,11,10,11,10,11,10,11];% these should also make an area with the same distances as the target sites
moveList_x = [];
moveList_y = [];
fakeList_x = [];
fakeList_y = [];
nspots = length(visCtx_x);
stim_rate = 1/40;%stimulation rate of 40Hz
stimProb= 0.5; %probability at which stimulation occurs
for i =1:length(visCtx_x)    
    moveList_x = [moveList_x,1, repmat(visCtx_x(i)*xConv,1,s2.Rate*stim_rate/nspots)];%generates lists for the galvos to move such that each of 8 spots is stimulatedat 40 hz
    moveList_y = [moveList_y,1, repmat(visCtx_y(i)*yConv,1,s2.Rate*stim_rate/nspots)];
    fakeList_x = [fakeList_x,1,repmat(offTarget_x(i)*xConv,1,s2.Rate*stim_rate/nspots)];
    fakeList_y = [fakeList_y,1,repmat(offTarget_y(i)*yConv,1,s2.Rate*stim_rate/nspots)];
end
moveList_x = repmat(moveList_x,1,300);
moveList_y = repmat(moveList_y,1,300);%pregenerated list of targets for stimulation on the mouse
fakeList_x = repmat(fakeList_x,1,300);
fakeList_y = repmat(fakeList_y,1,300);%pregenerated list that is off the mouse to use as false stimulation
%% Main execution loop    
XY_list = double.empty(0,length(offTarget_x)*2);

while true
if rand(1) <= stimProb %choose where the stimulation will be based on probability
    outputx = moveList_x';
    dest_x = visCtx_x;
    outputy = moveList_y';
    dest_y = visCtx_y;
else
    outputx = fakeList_x';
    dest_x = offTarget_x;
    outputy = fakeList_y';
    dest_y = offTarget_y;
end
    
  if ~ishandle(ButtonHandle)
      release(s0)
    disp('experiment stopped by user');
    break;
  end
  tic
queueOutputData(s2,[outputy outputx])%swap x and y because mirrors are swapped
newTrialListener = addlistener(s0,'DataAvailable', @newTrialCheckTargetted); %Add listener to check if there is a new trial aka if the laser should move 
s0.startBackground
toc
try
    s0.wait(61)%hold all operations for 61 seconds, (is more than max trial length), if new trial not triggered by then, this trial is repeated
    fprintf('Triggered!')
    XY_list = [XY_list; [dest_x,dest_y]];
catch %if 61 seconds elapses
    s2.startBackground 
    s0.stop()
    delete(newTrialListener); %delete(mirrorPosListener);
    disp("something went wrong, no trigger detected")
end
delete(newTrialListener); %delete(mirrorPosListener);

    lh = addlistener(s0,'DataAvailable',@stillRunningCheck);%if the scan is still running after a trial ends, stop it
    s0.startBackground
    
try
    s0.wait(7)
    s0.stop()
    delete(lh)
catch
    s0.stop()
    delete(lh)
    s2.stop()
end
end

%fclose(fid1);
%XY_list = XY_list(2:end,:);
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
    destListy = [repmat(3,6,1);repmat(2,6,1);repmat(1,8,1);repmat(0,8,1);repmat(-1,8,1);repmat(-2,8,1);repmat(-3,8,1);repmat(-4,8,1);repmat(-5,8,1);repmat(-6,4,1);repmat(-7,2,1)];
    destList = [destListx,destListy.*-1];%the y mirror is opposite so need to flip it
    idx=randsample([1:size(destList,1)],1,true);
    X = destList(idx,1); Y = destList(idx,2);
    if [X,Y] == XY_list(end,:)
        [X,Y]=getDestination(XY_list);
    end
end