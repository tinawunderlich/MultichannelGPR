function [traces,dt,ns,x,y,z]=readmala_single(foldername,name,use_GPS)

%
% Read Mala rd3-data
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
    temp=textscan(fid,'%f%s%s%f%s%f%s%f%s%f');
    fclose(fid);
    trnum=temp{1}; % Trace number
    for i=1:length(temp)
        if strcmp(temp{i}(1),'N')
            y=temp{i-1};
        end
        if strcmp(temp{i}(1),'E')
            x=temp{i-1};
        end
    end
    z=temp{8};
else
    fid=fopen(fullfile(foldername,[name,'.corc']),'r'); 
    temp=textscan(fid,'%f%s%s%f%s%f%s%f%s%f');
    fclose(fid);
    trnum=temp{1}; % Trace number
    for i=1:length(temp)
        if strcmp(temp{i}(1),'Y')
            y=temp{i-1};
        end
        if strcmp(temp{i}(1),'X')
            x=temp{i-1};
        end
    end
    z=temp{8};
end

%%% Read radar data file
fid=fopen(fullfile(foldername,[name,'.rd3']),'r');
traces=zeros(ns,num_traces);
for ii=1:num_traces    % for all traces in this file (also without coords)
    traces(:,ii)=fread(fid,ns,'int16','l');    % read trace
end
fclose(fid);

dt=range/ns;