function [traces,dt,ns,x,y,z,numchannels]=readmala_withoutGPS(foldername,name,profile_num)

% Read Mala rd7-data (each channel in one file)
% [traces,dt,ns,x,y,z,numchannels]=readmala_withoutGPS(foldername,name,profile_num)
%
% Dr. Tina Wunderlich, CAU Kiel 2024, tina.wunderlich@ifg.uni-kiel.de
%
% Input: 
% foldername: complete rSlicer-foldername and path
% name: Name of datafiles without '_...'
% profile_num: Numbers of profiles
% changeDir: flag if you want to change the +/- of the
% y-antenna-GPS-offset, yes=1, no=0
% add_Yoffset: additional offset for y-antenna-GPS-offset
%
% Output:
% traces is matrix of size [ns,numchannels*numtraces_per_channel] 
% (order is: first all traces of 1. channels, all traces of 2. channel,....)
% x, y, z: vectors of length [1,numchannels*numtraces_per_channel] with
% coordinates for each trace
% dt: sampling interval in ns
% ns: number of samples per trace
% numchannels: number of channels
%
% requires helmert.m


%%% Read yml-file
if profile_num<10
    temp=readlines(fullfile(foldername,[name,'_00',int2str(profile_num),'.yml']));
elseif profile_num>=10 && profile_num<100
    temp=readlines(fullfile(foldername,[name,'_0',int2str(profile_num),'.yml']));
else
    temp=readlines(fullfile(foldername,[name,'_',int2str(profile_num),'.yml']));
end
% get channel infos:
ind=find(strcmp(temp,'SELECTED CHANNELS:')); % find this line
anz=1;
while ~strcmp(temp(ind+anz),'TRANSMITTERS:')
    ch{anz}=extractAfter(temp(ind+anz),'- ');
    anz=anz+1;
end
test=cellfun(@(x) strsplit(x,':'),ch,'UniformOutput',false); % split R* and T*
channelinfo=zeros(length(test),3); % channel# T# R#
channelinfo(:,1)=1:length(test); % channel number
channelinfo(:,2)=cellfun(@(x) str2num(extractAfter(x{1},'T')),test); % transmitter #
channelinfo(:,3)=cellfun(@(x) str2num(extractAfter(x{2},'R')),test); % receiver #

% get positions of transmitters:
ind=find(strcmp(temp,'TRANSMITTERS:')); % find this line
i=1;
while ~strcmp(temp(ind+(i-1)*3+1),'RECIVERS:')
    Tx(i,1)=str2num(extractAfter(temp(ind+(i-1)*3+2),'X-POS: '));
    Ty(i,1)=str2num(extractAfter(temp(ind+(i-1)*3+3),'Y-POS: '));
    i=i+1;
end

% get positions of receivers:
ind=find(strcmp(temp,'RECIVERS:')); % find this line
i=1;
while ~strcmp(temp(ind+(i-1)*3+1),'WHEEL:')
    Rx(i,1)=str2num(extractAfter(temp(ind+(i-1)*3+2),'X-POS: '));
    Ry(i,1)=str2num(extractAfter(temp(ind+(i-1)*3+3),'Y-POS: '));
    i=i+1;
end



%%% Read header file
if profile_num<10
    list=dir(fullfile(foldername,[name,'_00',int2str(profile_num),'_*.rad'])); % list of all rad-files for this profile (one for each channel)
elseif profile_num>=10 && profile_num<100
    list=dir(fullfile(foldername,[name,'_0',int2str(profile_num),'_*.rad']));
else
    list=dir(fullfile(foldername,[name,'_',int2str(profile_num),'_*.rad']));
end
if isempty(list)
    traces=[];
    dt=[];
    ns=[];
    x=[];
    y=[];
    z=[];
    numchannels=[];
    return;
end

for j=1:length(list)
    fid=fopen(fullfile(foldername,list(j).name),'r');

    if fid==-1
        traces=[];
        dt=[];
        ns=[];
        x=[];
        y=[];
        z=[];
        numchannels=[];
        return;
    end

    i=1;
    while ~feof(fid)
        temp{i}=fgetl(fid);

        if strfind(temp{i},'SAMPLES:')
            ns(j)=str2double(temp{i}(length('SAMPLES:')+1:end));
        elseif strfind(temp{i},'TIMEWINDOW:')
            range(j)=str2double(temp{i}(length('TIMEWINDOW:')+1:end));
        elseif strfind(temp{i},'DISTANCE INTERVAL:')
            dx(j)=str2double(temp{i}(length('DISTANCE INTERVAL:')+1:end));
        elseif strfind(temp{i},'LAST TRACE:')
            num_traces(j)=str2double(temp{i}(length('LAST TRACE:')+1:end));
        elseif strfind(temp{i},'CHANNELS:')
            numchannels(j)=str2double(temp{i}(length('CHANNELS:')+1:end));
        elseif strfind(temp{i},'CHANNEL CONFIGURATION:')
            config{j}=temp{i}(length('CHANNEL CONFIGURATION:')+1:end);
            test=strsplit(config{j},':');
            T(j)=str2num(extractAfter(test{1},'T')); % transmitter number
            R(j)=str2num(extractAfter(test{2},'R')); % receiver number
        end
        i=i+1;
    end
    fclose(fid);
end


%%% Read radar data file (rd7 -> 32 bit int)
if profile_num<10
    list=dir(fullfile(foldername,[name,'_00',int2str(profile_num),'_*.rd7'])); % list of all channel files for this profile
elseif profile_num>=10 && profile_num<100
    list=dir(fullfile(foldername,[name,'_0',int2str(profile_num),'_*.rd7']));
else
    list=dir(fullfile(foldername,[name,'_',int2str(profile_num),'_*.rd7']));
end

traces1=cell(numchannels(1),1); % initialize cell array (one cell for each channel)
for i=1:numchannels(1)
    traces1{i}=zeros(ns(1),num_traces(1));
end

for j=1:length(list) % go through each channel for this profile

    fid=fopen(fullfile(foldername,list(j).name),'r');

    for ii=1:num_traces(j)    % for all traces in this file
        temp1=fread(fid,ns(j),'int32','l');    % read trace
        traces1{j}(:,ii)=temp1;
    end
    fclose(fid);
end

% sort in one matrix
traces=zeros(ns(1),num_traces(1)*numchannels(1));
for i=1:numchannels(1)
    traces(:,(i-1)*num_traces(1)+i-(i-1):i*num_traces(1))=traces1{i};
end

dt=range(1)/ns(1); % in ns
ns=ns(1); % number of samples

%%% Read pos-file
if profile_num<10
    fid=fopen(fullfile(foldername,[name,'_00',int2str(profile_num),'_G01.pos']),'r'); % list of all channel files for this profile
elseif profile_num>=10 && profile_num<100
    fid=fopen(fullfile(foldername,[name,'_0',int2str(profile_num),'_G01.pos']),'r');
else
    fid=fopen(fullfile(foldername,[name,'_',int2str(profile_num),'_G01.pos']),'r');
end
temp=textscan(fid,'%f%f%f%f','HeaderLines',1);
fclose(fid);

coords=[temp{1} temp{2} temp{3} temp{4}]; % trnum x y z -> two lines=start and end of profile

% interpolate for each trace
if num_traces(1)>coords(2,1)
    disp('  Number of data traces is larger than number of traces in coordinate file. Extrapolating linearly!')
end
xtemp=interp1(coords(:,1),coords(:,2),coords(1,1):num_traces(1),'linear','extrap');
ytemp=interp1(coords(:,1),coords(:,3),coords(1,1):num_traces(1),'linear','extrap');
ztemp=interp1(coords(:,1),coords(:,4),coords(1,1):num_traces(1),'linear','extrap');


%%%  apply Tx-Rx offsets
for i=1:numchannels(1)
    pos=zeros(num_traces(i),2); % initialize

    % offsets for this channel
    tcx=Tx(T(i));
    tcy=Ty(T(i));
    rcx=Rx(R(i));
    rcy=Ry(R(i));
    ch_x=mean([tcx rcx]); % midpoint in x direction
    ch_y=mean([tcy rcy]); % midpoint in y direction
    anz=15; % number of points for direction determination
    for ii=1:length(xtemp)-anz
        dist=sqrt((xtemp(ii)-xtemp(ii+anz))^2+(ytemp(ii)-ytemp(ii+anz))^2);
        pos(ii,:)=helmert([ch_x ch_y],[0 0; 0 dist],[xtemp(ii) ytemp(ii); xtemp(ii+anz) ytemp(ii+anz)]);
    end
    anz1=anz;
    % calculation for the last few traces
    anz=anz-1;
    for ii=length(xtemp)-anz1+1:length(xtemp)-1
        dist=sqrt((xtemp(ii)-xtemp(ii+anz))^2+(ytemp(ii)-ytemp(ii+anz))^2);
        pos(ii,:)=helmert([ch_x ch_y],[0 0; 0 dist],[xtemp(ii) ytemp(ii); xtemp(ii+anz) ytemp(ii+anz)]);
        anz=anz-1;
    end
    % extrapolate for last trace
    pos(ii+1,1)=interp1([1 2],[pos(ii-1,1) pos(ii,1)],3,'linear','extrap');
    pos(ii+1,2)=interp1([1 2],[pos(ii-1,2) pos(ii,2)],3,'linear','extrap');

    % set coordinates
    x((i-1)*num_traces(i)+i-(i-1):i*num_traces(i))=pos(:,1);
    y((i-1)*num_traces(i)+i-(i-1):i*num_traces(i))=pos(:,2);
    z((i-1)*num_traces(i)+i-(i-1):i*num_traces(i))=ztemp;
end


numchannels=numchannels(1);