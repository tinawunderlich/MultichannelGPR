function [traces,dt,ns,x,y,z]=readmala_single(foldername,name,use_GPS)

%
% Read Mala rd3-data (all channels in one file)
% [traces,dt,ns,x,y,z,numchannels]=readmala(foldername,name,use_GPS)
%
% Dr. Tina Wunderlich, CAU Kiel, June 2025, tina.wunderlich@ifg.uni-kiel.de
%
% requires helmert.m
%
% Input:
% foldername: path of folder
% name: Name of datafile
% use_GPS: =1 if with GPS, else =0
%
%
% Output:
% traces: matrix with traces of one profile
% dt: sampling interval in ns
% ns: Number of samples
% x,y,z: matrices with coordinates in m


%%% Read header file
if use_GPS==1
    fid=fopen(fullfile(foldername,[name,'.rad']),'r');
else
    fid=fopen(fullfile(foldername,[name,'.RAD']),'r');
end

if fid==-1
    traces=[];
    dt=[];
    ns=[];
    x=[];
    y=[];
    z=[];
    return;
end

i=1;
while ~feof(fid)
    temp{i}=fgetl(fid);
    if strfind(temp{i},'SAMPLES:')
        ns=str2double(temp{i}(length('SAMPLES:')+1:end));
    elseif strfind(temp{i},'TIMEWINDOW:')
        range=str2double(temp{i}(length('TIMEWINDOW:')+1:end));
    elseif strfind(temp{i},'LAST TRACE:')
        num_traces=str2double(temp{i}(length('LAST TRACE:')+1:end));
    end
    i=i+1;
end
fclose(fid);


%%% Coordinates
if use_GPS==1
    fid=fopen(fullfile(foldername,[name,'.cor']),'r');

fclose(fid);
else
    fid=fopen(fullfile(foldername,[name,'.corc']),'r');
    
fclose(fid);
end

%%% Read radar data file
if profile_num<10
    fid=fopen(fullfile(foldername,[name,'_00',int2str(profile_num),'.rd3']),'r');
elseif profile_num>=10 && profile_num<100
    fid=fopen(fullfile(foldername,[name,'_0',int2str(profile_num),'.rd3']),'r');
else
    fid=fopen(fullfile(foldername,[name,'_',int2str(profile_num),'.rd3']),'r');
end

traces=cell(numchannels,1); % initialize cell array (one cell for each channel)

for i=1:numchannels
    traces{i}=zeros(ns,num_traces/numchannels);
end
anz=1;
for ii=1:num_traces    % for all traces in this file (also without coords)
    temp1=fread(fid,ns,'int16','l');    % read trace
    currChan=mod(ii-1,numchannels)+1;  % current channel number
    traces{currChan}(:,anz)=temp1;
    if currChan==numchannels
        anz=anz+1;
    end
end
fclose(fid);

% delete traces without coords:
in=pos_new(:,1); % valid trace numbers with coordinates
for i=1:numchannels
    traces{i}=traces{i}(:,in);
end

dt=range/ns;